"""Import modules."""
from opentrons import protocol_api

metadata = {
    'protocolName': 'CPL protocol 10: Serial dilution. CPL10_230816_RUN1_TEMP_FIX',
    'author': 'Jacob Sturgess',
    'description': 'Serially dilute samples to a specified titre. (Uses 10ul into 90ul)',
    'apiLevel': '2.9'
}

def run(ctx):
    """Run the protocol."""
    # ------------------------------------------------------------------------------------------------------------------
    # USER DEFINED VARIABLES
    samples_to_dilute = 16
    dilutions_per_sample = 6 # Use a multiple of 12 to avoid remainder wells on rows
    include_neat = True # Should the first well in each dilution group be neat sample? (True/False)

    start_sm_buffer_volume = 16 # Volume to the nearest ml, for pipette height purposes
    # ------------------------------------------------------------------------------------------------------------------

    # Load Labware and pipettes
    dilution_plate = ctx.load_labware('nest_96_wellplate_200ul_flat', 1, 'Dilution plate')
    sample_plate = ctx.load_labware('nest_96_wellplate_200ul_flat', 2, 'Sample plate')
    tuberack = ctx.load_labware('opentrons_6_tuberack_falcon_50ml_conical', 3, 'Tube rack') # Position the rack so the tubes will be nearer the right hand side
    sm_buffer = tuberack.wells_by_name()['A1'] # The tube of SM buffer should go in the top left slot (A1)

    tiprack_300ul = ctx.load_labware('opentrons_96_tiprack_300ul', 4, '300ul Tip rack')
    p300 = ctx.load_instrument('p300_single_gen2', 'left', tip_racks=[tiprack_300ul])
    tiprack_20ul = ctx.load_labware('opentrons_96_tiprack_20ul', 5, '20ul Tip rack')
    p20 = ctx.load_instrument('p20_single_gen2', 'right', tip_racks=[tiprack_20ul])
    
    # Confirm user defined variables
    ctx.pause(f"""Please confirm these variables are correct:
                - Unique samples to dilute: {samples_to_dilute}.
                - Dilutions per sample: {dilutions_per_sample}.
                - Total wells produced: {samples_to_dilute * dilutions_per_sample}.
                - Start buffer volume: {start_sm_buffer_volume}.""")

    total_wells = samples_to_dilute * dilutions_per_sample
    sm_volume = start_sm_buffer_volume * 1000

    # Prepare plates for well-by-well iteration 
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

    dilution_wells = prepare_plate(dilution_plate, 96)
    sample_wells = prepare_plate(sample_plate, 96)

    dilution_plate_sublists = [dilution_wells[x:x+dilutions_per_sample] for x in range(0, total_wells, dilutions_per_sample)] # Create sublists in the dilution plate

    # Distribute SM buffer
    p300.pick_up_tip()
    if sm_volume > 50000: # Set pipette to aspirate just below the depth of the liquid in the 50ml tube
        p300.well_bottom_clearance.aspirate = 83 # This could all go in a for loop to be more accurate
    elif sm_volume > 45000:
        p300.well_bottom_clearance.aspirate = 74
    elif sm_volume > 40000:
        p300.well_bottom_clearance.aspirate = 65
    elif sm_volume > 35000:
        p300.well_bottom_clearance.aspirate = 56
    elif sm_volume > 30000:
        p300.well_bottom_clearance.aspirate = 47
    elif sm_volume > 25000:
        p300.well_bottom_clearance.aspirate = 38
    elif sm_volume > 20000:
        p300.well_bottom_clearance.aspirate = 29
    elif sm_volume > 15000:
        p300.well_bottom_clearance.aspirate = 20
    elif sm_volume > 10000:
        p300.well_bottom_clearance.aspirate = 11
    elif sm_volume > 5000:
        p300.well_bottom_clearance.aspirate = 3
    else:
        p300.well_bottom_clearance.aspirate = 1
    if include_neat == True:
        for sublist in dilution_plate_sublists:
            p300.distribute(90, source=sm_buffer, dest=sublist[1:], new_tip='never') # If dispensing neat sample, do not fill that well with SM buffer
    else:
        p300.distribute(90, source=sm_buffer, dest=dilution_wells[:total_wells], new_tip='never')
    p300.drop_tip()

    # Transfer neat samples
    for sample_well, dest_sublist in zip(sample_wells[:samples_to_dilute], dilution_plate_sublists):
        if include_neat == True:
            p20.transfer(20, source=sample_well, dest=dest_sublist[0]) # If we want to spot the neat sample, we will need some more liquid (e.g. 20ul in the well)
        else:
            p20.transfer(10, source=sample_well, dest=dest_sublist[0])

    # Dilute samples
    for dilution_sublist in dilution_plate_sublists:
        for i in range(len(dilution_sublist)):
            if dilution_sublist[i] != dilution_sublist[-1]:
                p20.transfer(10, source=dilution_sublist[i], dest=dilution_sublist[i+1], mix_after=(2, 20)) # Mix_after takes args: repetitions, volume

    ctx.comment(f'Protocol complete.')