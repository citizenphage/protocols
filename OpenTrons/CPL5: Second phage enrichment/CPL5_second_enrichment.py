"""Import modules."""
from opentrons import protocol_api

metadata = {
    'protocolName': 'CPL protocol 5: Second phage enrichment.',
    'author': 'Jacob Sturgess',
    'description': 'Adds 10ul host culture and 190ul growth media (LB + 0.25mM MgCl2 & CaCl2)to microplate. Then adds 5ul of filtered first enrichment. (Note: Each stage may be done independently.)',
    'apiLevel': '2.9'
}

def run(ctx):
    """Run the protocol."""
    # ------------------------------------------------------------------------------------------------------------------
    # USER DEFINED VARIABLES
    # Indicate how many wells of filtrate you are transferring and wells of reagents you require. Note: These steps can be done at different times if needed.
    filtrate_wells_to_transfer = 36
    wells_for_reagents = 36

    # Please also specify the rough volume (to the nearest ml) of lb and host culture that you are using for pipette height purposes:
    start_lb_volume = 20
    start_host_volume = 10
    # ------------------------------------------------------------------------------------------------------------------

    # Load labware and pipettes
    new_plate = ctx.load_labware('nest_96_wellplate_200ul_flat', 1, 'New microplate')
    if filtrate_wells_to_transfer > 0:
        first_filtrate = ctx.load_labware('nest_96_wellplate_200ul_flat', 2, '1st enrichment filtrate')
    if wells_for_reagents > 0:
        tuberack = ctx.load_labware('opentrons_6_tuberack_falcon_50ml_conical', 3, 'Tube rack') # Position the rack so the tubes will be nearer the right hand side
        lb = tuberack.wells_by_name()['A1'] # The tube of LB should go in the top left slot (A1)
        host = tuberack.wells_by_name()['B1'] # The tube of host culture should go in the bottom left slot (B1)
        p300tiprack = ctx.load_labware('opentrons_96_tiprack_300ul', 4, 'P300 Tip rack')
        p300 = ctx.load_instrument('p300_single_gen2', 'left', tip_racks=[p300tiprack])
    p20tiprack = ctx.load_labware('opentrons_96_tiprack_20ul', 5, 'P20 Tip rack')
    if filtrate_wells_to_transfer == 96:
        p20tiprack2 = ctx.load_labware('opentrons_96_tiprack_20ul', 6, 'P20 Tip rack 2')
        p20 = ctx.load_instrument('p20_single_gen2', 'right', tip_racks=[p20tiprack, p20tiprack2])
    else:
        p20 = ctx.load_instrument('p20_single_gen2', 'right', tip_racks=[p20tiprack])

    # Convert host volumes from ml to ul
    lb_volume, host_volume = start_lb_volume * 1000, start_host_volume * 1000

    # Checks
    if wells_for_reagents:
        if host_volume < wells_for_reagents*10 or lb_volume < wells_for_reagents*190:
            ctx.pause("Please ammend user defined variables. You may not have sufficient reagents.")

    # Confirm user defined variables
    ctx.pause(f"""Please confirm the user defined variables:
        filtrate_wells_to_transfer = {filtrate_wells_to_transfer}
        wells_for_reagents = {wells_for_reagents}
        start_lb_volume = {start_lb_volume}
        start_host_volume = {start_host_volume}
        If these are not correct, please exit protocol and ammend.""")
    
    # Enable iteration through microplates
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
    
    # Enable iteration through the empty microplate
    new_wells = prepare_plate(new_plate, 96)

    # Distribute reagents if inidicated by the user
    if wells_for_reagents > 0:
        ctx.comment("Distributing reagents...")
        p300.pick_up_tip()
        for well in range(0, wells_for_reagents, 2): # Distribute LB
            if lb_volume > 50000: # Set pipette to aspirate just below the depth of the liquid in the 50ml tube
                p300.well_bottom_clearance.aspirate = 86
            elif lb_volume > 45000:
                p300.well_bottom_clearance.aspirate = 77
            elif lb_volume > 40000:
                p300.well_bottom_clearance.aspirate = 68
            elif lb_volume > 35000:
                p300.well_bottom_clearance.aspirate = 59
            elif lb_volume > 30000:
                p300.well_bottom_clearance.aspirate = 50
            elif lb_volume > 25000:
                p300.well_bottom_clearance.aspirate = 41
            elif lb_volume > 20000:
                p300.well_bottom_clearance.aspirate = 32
            elif lb_volume > 15000:
                p300.well_bottom_clearance.aspirate = 23
            elif lb_volume > 10000:
                p300.well_bottom_clearance.aspirate = 14
            elif lb_volume > 5000:
                p300.well_bottom_clearance.aspirate = 5
            else:
                p300.well_bottom_clearance.aspirate = 1

            p300.distribute(150, source=lb, dest=new_wells[well:well+2], new_tip='never', carryover=False, disposal_volume=0)
            for well in new_wells[well:well+2]:
                lb_volume -= 150
        p300.drop_tip()

        p20.pick_up_tip()
        for well in range(0, wells_for_reagents, 2): # Distribute host culture
            if host_volume > 50000:
                p20.well_bottom_clearance.aspirate = 86
            elif host_volume > 45000:
                p20.well_bottom_clearance.aspirate = 77
            elif host_volume > 40000:
                p20.well_bottom_clearance.aspirate = 68
            elif host_volume > 35000:
                p20.well_bottom_clearance.aspirate = 59
            elif host_volume > 30000:
                p20.well_bottom_clearance.aspirate = 50
            elif host_volume > 25000:
                p20.well_bottom_clearance.aspirate = 41
            elif host_volume > 20000:
                p20.well_bottom_clearance.aspirate = 32
            elif host_volume > 15000:
                p20.well_bottom_clearance.aspirate = 23
            elif host_volume > 10000:
                p20.well_bottom_clearance.aspirate = 14
            elif host_volume > 5000:
                p20.well_bottom_clearance.aspirate = 5
            else:
                p20.well_bottom_clearance.aspirate = 1
            p20.distribute(10, source=host, dest=new_wells[well:well+2], new_tip='never', carryover=False, disposal_volume=0)
            for well in new_wells[well:well+2]:
                host_volume -= 10
        p20.drop_tip()
        ctx.comment(f"""
            Reagents were added to {wells_for_reagents} wells.
            LB volume went from {start_lb_volume*1000}ul to {lb_volume}ul.
            Host culture volume went from {start_host_volume*1000}ul to {host_volume}ul.""")

    # Transfer 5ul of each first filtrate to the fresh microplate, if indicated by the user
    if filtrate_wells_to_transfer > 0:
        filtrate_wells = prepare_plate(first_filtrate, 96)
        p20.well_bottom_clearance.aspirate = 1
        for filtrate_well, new_well in zip(filtrate_wells[:filtrate_wells_to_transfer], new_wells[:filtrate_wells_to_transfer]):
            p20.transfer(10, source=filtrate_well, dest=new_well)
            ctx.comment(f"{filtrate_wells_to_transfer} wells of first enrichment filtrate transfered.""")