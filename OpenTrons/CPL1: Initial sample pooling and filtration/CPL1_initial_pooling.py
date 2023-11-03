"""Import modules."""
from opentrons import protocol_api
import math

metadata = {
    'protocolName': 'CPL protocol 1: Initial pooling and filter sterilisation',
    'author': 'Jacob Sturgess',
    'description': 'Pool together wells of 100ul raw sample into a filter plate. For use with 1 or 2 raw wastewater plates. Can accomodate partially filled filter plate.',
    'apiLevel': '2.9'
}

def run(ctx):
    """Run the protocol."""
    # ------------------------------------------------------------------------------------------------------------------
    # USER DEFINED VARIABLES
    rows_in_rawplate1 = 8 # Number of rows with sample (1-8) in the 1st raw sample plate (deck position 2).
    rows_in_rawplate2 = 8 # Number of rows with sample (0-8) in the 2nd raw sample plate (deck position 3). Note - optional.
    full_rows_in_filterplate = 0 # Number of rows in the filter plate which already contain any liquid (0-8). Note - leave space to pool rawplate(s):
    # ------------------------------------------------------------------------------------------------------------------

    # Load Labware
    filterplate = ctx.load_labware('filterplate_96_wellplate_200ul', 1, 'Filter plate') # This is a custom labware definition
    rawplate1 = ctx.load_labware('nest_96_wellplate_200ul_flat', 2, 'Raw plate 1')
    if rows_in_rawplate2 != 0:
        rawplate2 = ctx.load_labware('nest_96_wellplate_200ul_flat', 3, 'Raw plate 2 (optional)')
    tiprack = ctx.load_labware('opentrons_96_tiprack_300ul', 4, 'Tip rack')
    # Load Pipette
    p300 = ctx.load_instrument('p300_single_gen2', 'right', tip_racks=[tiprack])

    # Check user variables
    if rows_in_rawplate1 > 8 or rows_in_rawplate1 <=0 or rows_in_rawplate2 > 8: 
        ctx.pause('Please check the user-defined values in protocol file.')
    if full_rows_in_filterplate + ((rows_in_rawplate1 + rows_in_rawplate2)/2) > 8:
        ctx.pause('Please check the values in protocol file. There may not be space for this many rows to be pooled')

    # Confirm user defined variables
    ctx.pause(f"""Please confirm the user defined variables:
                Rows_in_rawplate1 (deck position 2) = {rows_in_rawplate1}.
                Rows_in_rawplate2 (deck position 4) = {rows_in_rawplate2}.
                Full_rows_in_filterplate = {full_rows_in_filterplate}.
                If these are not correct, please exit protocol and ammend.""")

    def dispensing_raw(rawplate_row: int, dest_row: int, plate: int=1):
        """Dispense from the first rawplate.
        
        Args:
            rawplate_row: The row in the rawplate
            dest_row: Destination row in filterplate
            plate: Rawplate 1 or 2
        """
        if plate == 1:
            if rawplate_row >= rows_in_rawplate1: # Do not process if we have exceeded the available rows 
                return
            sources = rawplate1.rows()[rawplate_row]
            dests = filterplate.rows()[dest_row + full_rows_in_filterplate]

        elif plate == 2:
            if rawplate_row >= rows_in_rawplate2: # Do not process if we have exceeded the available rows 
                return
            sources = rawplate2.rows()[rawplate_row]
            dests = filterplate.rows()[dest_row + full_rows_in_filterplate + math.ceil(rows_in_rawplate1/2)]
        else:
            print("Incorrect plate number")
            return
        # Use the transfer function to move 100ul from the source wells to the destination wells
        p300.transfer(100, sources, dests)

    # Dispense rows of rawplate 1 into filterplate. Note: destination row indexes are repeated in pairs because 2 rows are combined. 
    ctx.comment('Pooling rawplate 1...')
    for src, dst in zip([0,1,2,3,4,5,6,7], [0,0,1,1,2,2,3,3]):
        dispensing_raw(src, dst, plate=1) 
    
    # Dispense rows of rawplate 2 into filterplate
    if rows_in_rawplate2 > 0:
        ctx.comment('Pooling rawplate 2...')
        for src, dst in zip([0,1,2,3,4,5,6,7], [0,0,1,1,2,2,3,3]):
            dispensing_raw(src, dst, plate=2) 

    ctx.comment('Pooling complete. Filter plate is ready to be sealed and sterilised by centifugation, then can be cold stored or used in "CPL2 Secondary Pooling" protocol.')