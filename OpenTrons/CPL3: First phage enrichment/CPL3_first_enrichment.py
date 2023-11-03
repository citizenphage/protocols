"""Import modules."""
from opentrons import protocol_api

metadata = {
    'protocolName': 'CPL protocol 3: First phage enrichment ',
    'author': 'Jacob Sturgess',
    'description': 'Adds 75ul host culture and 400ul growth media (4.4xLB + 0.25mM MgCl2 & CaCl2) with 1ml pools of filter-sterilised samples. Note: the samples may instead be added later.).',
    'apiLevel': '2.9'
}

def run(ctx):
    """Run the protocol."""
    # ------------------------------------------------------------------------------------------------------------------
    # USER DEFINED VARIABLES
    # Indicate how many wells you want to add growth media and host culture to:
    wells_for_reagents = 36

    # Please specify the rough volume (ml) of lb and host culture that you are using for pipette height purposes:
    start_lb_volume = 20
    start_host_volume = 20
    # ------------------------------------------------------------------------------------------------------------------

    # Load labware and pipettes
    deepwellplate = ctx.load_labware('nest_96_wellplate_2ml_deep', 1, 'Deep well plate')
    tuberack = ctx.load_labware('opentrons_6_tuberack_falcon_50ml_conical', 2, 'Tube rack') # Position the rack so the tubes will be nearer the right hand side
    lb = tuberack.wells_by_name()['A1'] # The tube of 4.4xLB should go in the top left slot (A1)
    host = tuberack.wells_by_name()['B1'] # The tube of host culture should go in the bottom left slot (B1)
    p300tiprack = ctx.load_labware('opentrons_96_tiprack_300ul', 4, 'P300 Tip rack')
    p1000tiprack = ctx.load_labware('opentrons_96_tiprack_1000ul', 5, 'P1000 Tip rack')
    p1000 = ctx.load_instrument('p1000_single_gen2', 'right', tip_racks=[p1000tiprack])
    p300 = ctx.load_instrument('p300_single_gen2', 'left', tip_racks=[p300tiprack])

    # Checks
    # Check sufficient volumes, number of wells...

    # Confirm user defined variables
    ctx.pause(f"""Please confirm the user defined variables:
        wells_for_reagents = {wells_for_reagents}
        start_lb_volume = {start_lb_volume}
        start_host_volume = {start_host_volume}
        If these are not correct, please exit protocol and ammend.""")

    # Enable iteration through deepwell plate
    deeps = []
    for row in range(8):
        for well in range(12):
            deeps.append(deepwellplate.rows()[row][well])

    # Set pipettes to dispense near the top of the deep well (35mm from base) to avoid liquid splashing the tip
    p1000.well_bottom_clearance.dispense, p300.well_bottom_clearance.dispense = 35, 35
    # Convert host volumes from ml to ul
    lb_volume, host_volume = start_lb_volume * 1000, start_host_volume * 1000 

    # Set source and destinations for distributing LB
    p1000.pick_up_tip()
    for well in range(0, wells_for_reagents, 2):
        source = lb
        dests = list(deeps[well:well+2])

        # Set pipette to aspirate just below the depth of the liquid in the 50ml tube
        if lb_volume > 50000:
            p1000.well_bottom_clearance.aspirate = 91
        elif lb_volume > 45000:
            p1000.well_bottom_clearance.aspirate = 82
        elif lb_volume > 40000:
            p1000.well_bottom_clearance.aspirate = 73
        elif lb_volume > 35000:
            p1000.well_bottom_clearance.aspirate = 64
        elif lb_volume > 30000:
            p1000.well_bottom_clearance.aspirate = 55
        elif lb_volume > 25000:
            p1000.well_bottom_clearance.aspirate = 46
        elif lb_volume > 20000:
            p1000.well_bottom_clearance.aspirate = 37
        elif lb_volume > 15000:
            p1000.well_bottom_clearance.aspirate = 28
        elif lb_volume > 10000:
            p1000.well_bottom_clearance.aspirate = 19
        elif lb_volume > 5000:
            p1000.well_bottom_clearance.aspirate = 10
        else:
            p1000.well_bottom_clearance.aspirate = 1.5

        # Distribute and update volume
        p1000.distribute(400, source, dests, new_tip='never', carryover = False, disposal_volume=0)
        for well in dests:
            lb_volume -= 400
    p1000.drop_tip()

    # Distribute 75ul of host to indicated wells
    p300.pick_up_tip()
    for well in range(0, wells_for_reagents, 4):
        source = host
        dests = list(deeps[well:well+4])

        # Set pipette to aspirate just below the depth of the liquid in the 50ml tube
        if host_volume > 50000:
            p300.well_bottom_clearance.aspirate = 91
        elif host_volume > 45000:
            p300.well_bottom_clearance.aspirate = 82
        elif host_volume > 40000:
            p300.well_bottom_clearance.aspirate = 73
        elif host_volume > 35000:
            p300.well_bottom_clearance.aspirate = 64
        elif host_volume > 30000:
            p300.well_bottom_clearance.aspirate = 55
        elif host_volume > 25000:
            p300.well_bottom_clearance.aspirate = 46
        elif host_volume > 20000:
            p300.well_bottom_clearance.aspirate = 37
        elif host_volume > 15000:
            p300.well_bottom_clearance.aspirate = 28
        elif host_volume > 10000:
            p300.well_bottom_clearance.aspirate = 19
        elif host_volume > 5000:
            p300.well_bottom_clearance.aspirate = 10
        else:
            p300.well_bottom_clearance.aspirate = 1.5

        # Distribute and update volume
        p300.distribute(75, source, dests, new_tip='never', carryover=False, disposal_volume=0)
        for well in dests:
            host_volume -= 75
    p300.drop_tip()

    ctx.comment(f"""
            Complete. Reagents were added to {wells_for_reagents} wells.
            LB volume went from {start_lb_volume*1000}ul to {lb_volume}ul.
            Host culture volume went from {start_host_volume*1000}ul to {host_volume}ul.
            Deep well plate ready for shaking incubator overnight.""")