"""Import modules."""
from opentrons import protocol_api
import math

metadata = {
    'protocolName': 'CPL protocol 2: Secondary Pooling',
    'author': 'Jacob Sturgess',
    'description': 'Pools 5 wells of 200ul sterilised sample into wells in the deep well plate. This plate can then be cold stored or taken through to first round of phage enrichment.',
    'apiLevel': '2.9'
}

def run(ctx):
    """Run the protocol."""
    # ------------------------------------------------------------------------------------------------------------------
    # USER DEFINED VARIABLES
    # Indicate a plate is not in use by setting it equal to 0. Load plates in the order below.

    sampleplate1 = 96       # Wells of sample in sample plate 1 (deck position 2)
    sampleplate2 = 0        # Wells of sample in sample plate 2 (deck position 3)
    sampleplate3 = 0        # Wells of sample in sample plate 3 (deck position 4)
    sampleplate4 = 0        # Wells of sample in sample plate 4 (deck position 5)
    sampleplate5 = 0        # Wells of sample in sample plate 5 (deck position 6)

    full_wells_deepwell = 0 # Wells already full in the deep well plate
    # ------------------------------------------------------------------------------------------------------------------

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
    
    # Initialise variables for destination well and pools size.
    start_well = 0
    pool_size = 5

    # Load Labware and pipette
    deepwellplate = ctx.load_labware('nest_96_wellplate_2ml_deep', 1, 'Deep well plate')
    deepwells = prepare_plate(deepwellplate, 96)
    deepwells = iter(deepwells[full_wells_deepwell:]) # Remove unavailable wells which already contain liquid.
    p1000tiprack = ctx.load_labware('opentrons_96_tiprack_1000ul', 7, ' p1000 Tip rack')
    p1000 = ctx.load_instrument('p1000_single_gen2', 'right', tip_racks = [p1000tiprack])

    # Only load plates with liquid and add them to all_samples.
    all_samples = []
    for idx, wells_in_plate in enumerate([sampleplate1, sampleplate2, sampleplate3, sampleplate4, sampleplate5]):
        if wells_in_plate != 0:
            samples = ctx.load_labware('nest_96_wellplate_200ul_flat', idx+2, f'Sample plate {idx+1}')
            all_samples += prepare_plate(samples, wells_in_plate)

    # Confirm user defined variables
    ctx.pause(f"""Please confirm the user defined variables:
    sampleplate1 = {sampleplate1}
    sampleplate2 = {sampleplate2}
    sampleplate3 = {sampleplate3}
    sampleplate4 = {sampleplate4}
    sampleplate5 = {sampleplate5}
    full_wells_deepwell = {full_wells_deepwell}
    If these are not correct, please exit protocol and ammend.""")

    # Consolidate every 5 samples (each of 200ul volume) from the sampleplates into the deep well plate.
    while start_well < len(all_samples):
        sources = all_samples[start_well : start_well + pool_size]
        dest = next(deepwells)
        p1000.consolidate(200, sources, dest)
        start_well += pool_size

    ctx.comment(f'Complete. {math.ceil(all_samples / 5)} pools created. Deep well ready for cold storage or first phage enrichment.')