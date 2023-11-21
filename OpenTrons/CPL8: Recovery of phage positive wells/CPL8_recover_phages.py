"""Import modules."""
from opentrons import protocol_api

metadata = {
    'protocolName': 'CPL protocol 8: Recovery of phage positive wells',
    'author': 'Jacob Sturgess',
    'description': 'Parse through the csv file produced by the “CPL R” pipeline, to identify the wells in the 2nd enrichment filtrate which will contain phages. Transfer the content of these wells to another microplate for downstream processing or other assays.',
    'apiLevel': '2.9'
}

def run(ctx):
    """Run the protocol."""
    # ------------------------------------------------------------------------------------------------------------------
    # USER DEFINED VARIABLES
    # Define the path to the csv produced by the "CPL R" pipeline
    phage_wells_csv = r"path/to/csv"
    # ------------------------------------------------------------------------------------------------------------------

    # Load labware and pipette
    second_filtrate = ctx.load_labware('nest_96_wellplate_200ul_flat', 1, '2nd enrichment filtrate')
    empty_microplate = ctx.load_labware('nest_96_wellplate_200ul_flat', 2, 'Microplate for phages')
    p300tiprack = ctx.load_labware('opentrons_96_tiprack_300ul', 4, 'P300 Tip rack')
    p300 = ctx.load_instrument('p300_single_gen2', 'left', tip_racks=[p300tiprack])

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
    
    microplate = prepare_plate(empty_microplate, 96)
    
    # Confirm contents of csv with the user
    with open(phage_wells_csv, 'r') as file: # Read csv file using a context manager
        ctx.comment("""Please check csv file details (showing first three rows): """)
        row_counter = 0
        for line in file:
            if row_counter == 3: # Set limit on rows displayed in preview
                break
            row = line.strip().split(',')
            ctx.comment(str(row))
            row_counter += 1
        ctx.pause("Confirm the above details before you continue.")

    with open(phage_wells_csv, 'r') as file: # Re-open csv file
        row_counter = 0
        for line, well in zip(file, microplate): # Loop through the csv and the microplate
            row = line.strip().split(',') # Remove leading/trailing whitespace and split by comma 
            if str(row[1]) == "Vp": # Skip the header row iteration
                continue
            # Set source and destination wells, and transfer wells indicated to contain phages
            source = second_filtrate[str(row[0])]
            dest = well
            p300.transfer(200, source, dest)
            row_counter += 1
        
    ctx.comment(f"Protocol complete. The microplate is ready to be sealed. {row_counter} wells have been transfered.")
