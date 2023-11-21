"""Import modules."""
from opentrons import protocol_api

metadata = {
    'protocolName': 'CPL protocol 9: Omnitray spot assay . CPL9_230816_RUN1',
    'author': 'Jacob Sturgess',
    'description': 'Spot wells in a 96 well plate onto an omnitray (which should have 45ml bottom agar as a base layer, with a mixture of 4.5ml top agar + 1.5ml host poured evenly on top).',
    'apiLevel': '2.9'
}

def run(ctx):
    """Run the protocol."""
    # ------------------------------------------------------------------------------------------------------------------
    # USER DEFINED VARIABLES
    # Classic spot assay
    wells_to_spot = 36
    use_alternate_spotting = True # Space out spots checkerboard-style, only possible with half a plate of samples or less (True/False)

    # Extra configurations
    volume_to_spot = 4 # Default = 4 (ul)
    dispense_height = 4 # Default = 4 (mm) when using 40ml of bottom agar + 4.5ml top agar + 1.5ml host. Refers to the height above the omnitray well.
    test_height_before_starting = True # Run a confirmation of the dispense height at protocol start (True/False)

    # Optional serial dilution configuration
    serial_dilution_spotting = False # Default = False. Can help save tips by moving up in titre within dilutions of the same sample
    num_dilutions_per_sample = 6 # Use either 3, 4, 6 or 12 dilutions per sample to fit on a row without remainder
    num_samples_used = 16 # For example; maximum = 16 when 6 dilutions per sample used
    # ------------------------------------------------------------------------------------------------------------------

    # Load Labware and pipette
    omnitray_96 = ctx.load_labware('nunc.omnitray_96_wellplate_20ul', 1, 'Omnitray') # Custom omnitray definition designed for spotting 96 wells onto
    well_plate = ctx.load_labware('nest_96_wellplate_200ul_flat', 2, '96 well plate')
    tiprack_20ul = ctx.load_labware('opentrons_96_tiprack_20ul', 4, '20ul Tip rack')
    p20 = ctx.load_instrument('p20_single_gen2', 'right', tip_racks=[tiprack_20ul])

    # Checks
    if use_alternate_spotting:
        if wells_to_spot > 48:
            ctx.pause("Cannot use alternate spotting for this many wells.")
    if serial_dilution_spotting:
        if num_dilutions_per_sample * num_samples_used > 96 or num_dilutions_per_sample * num_samples_used <= 0:
            ctx.pause("Amend the variables for serial dilution spotting.")

    # Confirm user defined variables
    ctx.pause(f"""Please confirm these variables are correct:
                wells_to_spot = {wells_to_spot}.
                use_alternate_spotting = {use_alternate_spotting}.

                volume_to_spot = {volume_to_spot}.
                dispense_height = {dispense_height}.
                test_height_before_starting = {dispense_height}.

                serial_dilution_spotting = {serial_dilution_spotting}.
                num_dilutions_per_sample = {num_dilutions_per_sample}.
                num_samples_used = {num_samples_used}.""")

    # Dispense height check
    if test_height_before_starting:
        p20.pick_up_tip()
        p20.move_to(omnitray_96['A1'].bottom(z=dispense_height))
        ctx.pause(f"Confirm this height ({dispense_height}) is correct.")
        p20.drop_tip()

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
 
    omnitray_spot_locations = prepare_plate(omnitray_96, 96)
    plate_wells = prepare_plate(well_plate, 96)

    p20.well_bottom_clearance.dispense = dispense_height

    # Classic spotting
    if serial_dilution_spotting == False and wells_to_spot != 0:

        if use_alternate_spotting == True and wells_to_spot <= 48: # Alternate spotting only possible with 48 or less wells
            row_counter = 1
            omnitray_spot_locations = [] # Overwrite omnitray spot locations list
            for row in omnitray_96.rows(): # For every row in the plate
                if row_counter % 2 != 0: # If it is an odd row (1,3,5,7,9,11)
                    for num in range(0,11,2): 
                        omnitray_spot_locations.append(row[num]) # Add every other well, starting at the first.
                else: # If it is an even row (2,4,6,8,10,12)
                    for num in range(1,12,2):
                        omnitray_spot_locations.append(row[num]) # Add every other well, starting at the second. 
                row_counter += 1

        for location, well in zip(omnitray_spot_locations[:wells_to_spot], plate_wells[:wells_to_spot]):
            p20.transfer(volume_to_spot, source=well, dest=location, blow_out=True, blowout_location='destination well')
        
        ctx.comment(f'Protocol complete. {wells_to_spot} wells were spotted onto the omnitray.')


    # Serial dilution spotting 
    if serial_dilution_spotting == True and num_dilutions_per_sample != 0 and num_samples_used != 0:

        well_total = num_dilutions_per_sample * num_samples_used
        samples = [plate_wells[x:x+num_dilutions_per_sample] for x in range(0, well_total, num_dilutions_per_sample)] # Create sublists of samples
        locations = [omnitray_spot_locations[y:y+num_dilutions_per_sample] for y in range(0, well_total, num_dilutions_per_sample)] # Create sublists of omnitray locations to match

        for sample_sublist, loc_sublist in zip(samples, locations):
            sample_sublist.reverse() # Reverse sublists so samples are spotted from low titre to high, so same tip can be used
            loc_sublist.reverse()
            p20.pick_up_tip()
            for sample, loc in zip(sample_sublist, loc_sublist):
                p20.transfer(volume_to_spot, source=sample, dest=loc, new_tip='never', blow_out=True, blowout_location='destination well')
            p20.drop_tip()
        
        ctx.comment(f'Protocol complete. {num_samples_used * num_samples_used} wells were spotted onto the omnitray.')
