
from opentrons import protocol_api
import math

metadata = {
    'protocolName': 'Sample pooling',
    'author': 'Jacob Sturgess',
    'description': 'Pools together wells of 100ul raw sample into a filter plate. Can accomodate 1 or 2 raw plates, and partially filled filter plate. User defined variables are "rows_in_rawplate1", "rows_in_rawplate2" and "full_rows_in_filterplate."',
    'apiLevel': '2.9'
}

def run(ctx):


    # ------------------------------------------------------------------------------------------------------------------
    # USER DEFINED VARIABLES
    # Number of rows with sample (1-8) in the 1st raw sample plate (at deck position 2):
    rows_in_rawplate1 = 3
    # Number of rows with sample (0-8) in the 2nd raw sample plate (at deck position 4). Note - this plate is optional:
    rows_in_rawplate2 = 1
    # Number of rows in the filter plate which already contain any liquid (0-8). Note - leave space to pool rawplate(s):
    full_rows_in_filterplate = 0
    # ------------------------------------------------------------------------------------------------------------------




    # Load Labware
    filterplate = ctx.load_labware('filterplate_96_wellplate_200ul', 1, 'Filter plate') # This is a custom labware definition
    rawplate1 = ctx.load_labware('nest_96_wellplate_200ul_flat', 2, 'Raw plate 1')
    rawplate2 = ctx.load_labware('nest_96_wellplate_200ul_flat', 4, 'Raw plate 2 (optional)') # Optional raw plate
    tiprack = ctx.load_labware('opentrons_96_tiprack_300ul', 3, 'Tip rack')

    # Load Pipette
    #m300 = ctx.load_instrument('p300_multi_gen2', 'left', tip_racks=[tiprack])
    p300 = ctx.load_instrument('p300_single_gen2', 'right', tip_racks=[tiprack])


    # Checks
    if rows_in_rawplate1 > 8: 
        ctx.pause(f'Please check the "rows_in_rawplate1" value in protocol file. It currently reads {rows_in_rawplate1}.')
    if rows_in_rawplate1 <=0: 
        ctx.pause(f'Please check the "rows_in_rawplate1" value in protocol file. It currently reads {rows_in_rawplate1}.')
    if rows_in_rawplate2 > 8: 
        ctx.pause(f'Please check the "rows_in_rawplate2" value in protocol file. It currently reads {rows_in_rawplate2}.')
    if rows_in_rawplate1 + rows_in_rawplate2 > 16:
        ctx.pause(f'Please check the "rows_in_rawplate1" and "rows_in_rawplate2" values in protocol file. There is not space for this many rows of sample to be pooled')
    if full_rows_in_filterplate + ((rows_in_rawplate1 + rows_in_rawplate2)/2) > 8:
        ctx.pause(f'Please check the values in protocol file. Ensure there is space for this many rows of sample to be pooled')

    # Confirm user defined variables
    ctx.pause(f'Please confirm there are {rows_in_rawplate1} rows of sample in the near raw plate (deck position 2) and {rows_in_rawplate2} rows in the far raw plate (deck position 4). Please also confirm there are {full_rows_in_filterplate} rows already containing liquid in the filter plate. If these are not correct, please exit protocol and edit the user defined variables.')


    ctx.comment('Pooling rawplate 1...')
    # Move row A (1st) of rawplate (100ul) into filterplate. Row indexing considers already full rows in filter plate.

    if rows_in_rawplate1 >= 1:
        sources = rawplate1.rows()[0]
        newindex = 0 + full_rows_in_filterplate
        dests = filterplate.rows()[newindex]
        p300.transfer(100, sources, dests)

    # Row B (2nd) of rawplate into filterplate
    if rows_in_rawplate1 >= 2:
        sources = rawplate1.rows()[1]
        newindex = 0 + full_rows_in_filterplate
        dests = filterplate.rows()[newindex]
        p300.transfer(100, sources, dests)

    # Row C (3rd) of rawplate into filterplate
    if rows_in_rawplate1 >= 3:
        sources = rawplate1.rows()[2]
        newindex = 1 + full_rows_in_filterplate
        dests = filterplate.rows()[newindex]
        p300.transfer(100, sources, dests)

    # Row D (4th) of rawplate into filterplate
    if rows_in_rawplate1 >= 4:
        sources = rawplate1.rows()[3]
        newindex = 1 + full_rows_in_filterplate
        dests = filterplate.rows()[newindex]
        p300.transfer(100, sources, dests)

    # Row E (5th) of rawplate into filterplate
    if rows_in_rawplate1 >= 5:
        sources = rawplate1.rows()[4]
        newindex = 2 + full_rows_in_filterplate
        dests = filterplate.rows()[newindex]
        p300.transfer(100, sources, dests)

    # Row F (6th) of rawplate into filterplate
    if rows_in_rawplate1 >= 6:
        sources = rawplate1.rows()[5]
        newindex = 2 + full_rows_in_filterplate
        dests = filterplate.rows()[newindex]
        p300.transfer(100, sources, dests)

    # Row G (7th) of rawplate into filterplate
    if rows_in_rawplate1 >= 7:
        sources = rawplate1.rows()[6]
        newindex = 3 + full_rows_in_filterplate
        dests = filterplate.rows()[newindex]
        p300.transfer(100, sources, dests)

    # Row H (8th) of rawplate into filterplate
    if rows_in_rawplate1 >= 8:
        sources = rawplate1.rows()[7]
        newindex = 3 + full_rows_in_filterplate
        dests = filterplate.rows()[newindex]
        p300.transfer(100, sources, dests)
    
    # Pooling actions for rawplate 2 (optional plate)
    if rows_in_rawplate2 > 0:
        ctx.comment('Pooling rawplate 2...')

    # Move row A (1st) of rawplate2 (100ul) into filterplate. Row indexing considers already full rows in filter plate and pooled rows from rawplate1.
    # Using floor() and ceil() to ensure rows that only contain 100ul of rawplate1 samples can be properly filled.
    if rows_in_rawplate2 >= 1:
        sources = rawplate2.rows()[0]
        newindex = 0 + full_rows_in_filterplate + math.floor(rows_in_rawplate1/2)
        dests = filterplate.rows()[newindex]
        p300.transfer(100, sources, dests)

    # Row B (2nd) of rawplate2 into filterplate
    if rows_in_rawplate2 >= 2:
        sources = rawplate2.rows()[1]
        newindex = 0 + full_rows_in_filterplate + math.ceil(rows_in_rawplate1/2)
        dests = filterplate.rows()[newindex]
        p300.transfer(100, sources, dests)

    # Row C (3rd) of rawplate2 into filterplate
    if rows_in_rawplate2 >= 3:
        sources = rawplate2.rows()[2]
        newindex = 1 + full_rows_in_filterplate + math.floor(rows_in_rawplate1/2)
        dests = filterplate.rows()[newindex]
        p300.transfer(100, sources, dests)

    # Row D (4th) of rawplate2 into filterplate
    if rows_in_rawplate2 >= 4:
        sources = rawplate2.rows()[3]
        newindex = 1 + full_rows_in_filterplate + math.ceil(rows_in_rawplate1/2)
        dests = filterplate.rows()[newindex]
        p300.transfer(100, sources, dests)

    # Row E (5th) of rawplate2 into filterplate
    if rows_in_rawplate2 >= 5:
        sources = rawplate2.rows()[4]
        newindex = 2 + full_rows_in_filterplate + math.floor(rows_in_rawplate1/2)
        dests = filterplate.rows()[newindex]
        p300.transfer(100, sources, dests)

    # Row F (6th) of rawplate2 into filterplate
    if rows_in_rawplate2 >= 6:
        sources = rawplate2.rows()[5]
        newindex = 2 + full_rows_in_filterplate + math.ceil(rows_in_rawplate1/2)
        dests = filterplate.rows()[newindex]
        p300.transfer(100, sources, dests)

    # Row G (7th) of rawplate2 into filterplate
    if rows_in_rawplate2 >= 7:
        sources = rawplate2.rows()[6]
        newindex = 3 + full_rows_in_filterplate + math.floor(rows_in_rawplate1/2)
        dests = filterplate.rows()[newindex]
        p300.transfer(100, sources, dests)

    # Row H (8th) of rawplate2 into filterplate
    if rows_in_rawplate2 >= 8:
        sources = rawplate2.rows()[7]
        newindex = 3 + full_rows_in_filterplate + math.ceil(rows_in_rawplate1/2)
        dests = filterplate.rows()[newindex]
        p300.transfer(100, sources, dests)


    ctx.comment('Pooling complete. Filter plate is ready to be sealed and sterilised by centifugation, then can be cold stored or used in "First phage enrichment" protocol.')
    
