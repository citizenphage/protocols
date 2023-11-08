"""Import modules."""
from opentrons import protocol_api

metadata = {
    'protocolName': 'CPL protocol 7: Microplate assay.',
    'author': 'Jacob Sturgess',
    'description': 'Adds LB, host and 2nd enrichment filtrate to microplate to run in plate reader assay - measuring optical density in each well over time. Includes addition of controls.',
    'apiLevel': '2.9'
}

def run(ctx):
    """Run the protocol."""
    # ------------------------------------------------------------------------------------------------------------------
    # USER DEFINED VARIABLES
    filtrate_wells_to_transfer = 36 # Indicate how many wells of filtrate you are transferring. Note: Maximum is 94 becasuse the assay requires space for 2 control wells

    # Please also specify the rough volume (to the nearest ml) of lb and host culture  that you are using for pipette height purposes:
    start_lb_volume = 30
    start_host_volume = 30
    # ------------------------------------------------------------------------------------------------------------------

    # Load labware and pipettes
    microplate = ctx.load_labware('nest_96_wellplate_200ul_flat', 1, 'New microplate')
    second_filtrate = ctx.load_labware('nest_96_wellplate_200ul_flat', 2, '2nd enrichment filtrate')
    tuberack = ctx.load_labware('opentrons_6_tuberack_falcon_50ml_conical', 3, 'Tube rack') # Position the rack so the tubes will be nearer the right hand side
    lb = tuberack.wells_by_name()['A1'] # The tube of LB should go in the top left slot (A1)
    host = tuberack.wells_by_name()['B1'] # The tube of host culture should go in the top centre slot (A2)
    p300tiprack = ctx.load_labware('opentrons_96_tiprack_300ul', 4, 'P300 Tip rack')
    p300 = ctx.load_instrument('p300_single_gen2', 'left', tip_racks=[p300tiprack])
    p20tiprack = ctx.load_labware('opentrons_96_tiprack_20ul', 5, 'P20 Tip rack')
    p20 = ctx.load_instrument('p20_single_gen2', 'right', tip_racks=[p20tiprack])

    neg_ctrl_well = microplate['H11']
    blank_well = microplate['H12']

    # Checks
    required_lb_ul = (filtrate_wells_to_transfer * 90) + 200
    required_host_ul = filtrate_wells_to_transfer * 5
    if start_host_volume*1000 < required_host_ul or start_lb_volume*1000 < required_lb_ul:
        ctx.pause("You may require larger volumes of LB or host to run this assay.")
    if filtrate_wells_to_transfer > 94:
        ctx.pause("Protocol limit for filtrate wells is 94.")

    # Confirm user defined variables
    ctx.pause(f"""Please confirm these variables are accurate:
            Wells of filtrate: {filtrate_wells_to_transfer},
            LB volume: {start_lb_volume},
            Host culture volume: {start_host_volume}.""")

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
    
    # Prepare plates for iteration
    new_wells = prepare_plate(sample=microplate, num_wells=filtrate_wells_to_transfer) 
    filtrate_wells = prepare_plate(sample=second_filtrate, num_wells=filtrate_wells_to_transfer)

    # Convert host volumes from ml to ul
    lb_volume, host_volume = start_lb_volume * 1000, start_host_volume * 1000

    # Distribute LB
    p300.pick_up_tip()
    for well in range(0, filtrate_wells_to_transfer, 3):
        if lb_volume > 50000:
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

        p300.distribute(90, source=lb, dest=new_wells[well:well+3], new_tip='never', carryover=False, disposal_volume=0)
        for well in new_wells[well:well+3]:
            lb_volume -= 90
    p300.transfer(95, source=lb, dest=neg_ctrl_well, new_tip='never') # Add lb to control wells
    p300.transfer(100, source=lb, dest=blank_well, new_tip='never')
    lb_volume -= 195
    p300.drop_tip()

    # Distribute host 
    p20.pick_up_tip()
    for well in range(0, filtrate_wells_to_transfer, 4): 
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

        p20.distribute(5, source=host, dest=new_wells[well:well+4], new_tip='never', carryover=False, disposal_volume=0)
        for well in new_wells[well:well+4]:
            host_volume -= 5
    p20.transfer(5, source=host, dest=neg_ctrl_well, new_tip='never') # Add host to negative control well
    host_volume -= 5
    p20.drop_tip()

    # Transfer 5ul of the second enrichment filtrate to the microplate
    p20.well_bottom_clearance.aspirate = 1
    for filtrate_well, new_well in zip(filtrate_wells, new_wells):
        p20.transfer(5, source=filtrate_well, dest=new_well)

    ctx.comment(f"""
        Complete.
        {filtrate_wells_to_transfer} wells of filtrate were added for the microplate assay.
        LB volume went from {start_lb_volume*1000}ul to {lb_volume}ul.
        Host culture volume went from {start_host_volume*1000}ul to {host_volume}ul.
        Well 'H11' is the negative control well. Well 'H12' is the blank well.""")
