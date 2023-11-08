"""Import modules."""
from opentrons import protocol_api

metadata = {
    'protocolName': 'CPL protocol 6: Filter second enrichment.',
    'author': 'Jacob Sturgess',
    'description': 'Transfer 200ul from each well of second enrichment to a filter plate.',
    'apiLevel': '2.9'
}

def run(ctx):
    """Run the protocol."""
    # ------------------------------------------------------------------------------------------------------------------
    # USER DEFINED VARIABLES
    wells_to_filter = 36
    full_wells_in_filterplate = 0 # Number of wells in the filter plate which already contain any liquid.
    # ------------------------------------------------------------------------------------------------------------------

    # Load Labware and pipette
    enrichment_plate = ctx.load_labware('nest_96_wellplate_200ul_flat', 1, '2nd enrichment plate')
    filterplate = ctx.load_labware('filterplate_96_wellplate_200ul', 2, 'Filter plate') # Custom labware definition
    tiprack = ctx.load_labware('opentrons_96_tiprack_300ul', 4, 'Tip rack')
    p300 = ctx.load_instrument('p300_single_gen2', 'left', tip_racks=[tiprack])

    # Checks
    if wells_to_filter <= 0 or wells_to_filter > 96 or wells_to_filter + full_wells_in_filterplate > 96:
        ctx.pause("Please ammend user defined variables.")

    # Confirm user defined variables
    ctx.pause(f"""Please confirm the user defined variables:
        wells_to_filter = {wells_to_filter}
        full_wells_in_filterplate = {full_wells_in_filterplate}
        If these are not correct, please exit protocol and ammend.""")

    # Main pipetting actions
    def prepare_plate(sample, num_wells: int):
        """Flatten a 2D array of wells into a 1d array.
        
        Args:
            sample: The labware form of the sampleplate.
            num_wells: User defined quantity of wells.
        """
        wells = []
        for row in sample.rows():
            if num_wells < 12:
                for well in row[:num_wells]:
                    wells.append(well)
                break
            for well in row:
                wells.append(well)
            num_wells -= 12
        return wells

    enrichment_wells = prepare_plate(enrichment_plate, 96)[:wells_to_filter]
    filterwells = prepare_plate(filterplate, 96)[full_wells_in_filterplate : wells_to_filter]

    for well, filterwell in zip(enrichment_wells, filterwells):
        p300.transfer(200, well, filterwell)

    ctx.comment('Protocol complete. Filter plate is ready to be sealed and sterilised by centifugation, then can be cold stored or used in the "CPL7 microplate assay" protocol.')