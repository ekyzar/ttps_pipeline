#**BrainJ Single Brain Processing**
#
#- Currently expects registered sections from BrainJ ImageJ pipeline
#- Channels are enhanced using CARE networks trained on 1.6µm data and segmented using U-Nets and StarDist.
#- csv files of raw, transformed and mapped cells are generated in addition to optional intermediate images for data checking

#Add the BrainJ install directory to the path
import sys
import logging
from datetime import datetime
sys.path.append('C:/Users/Alejandro/Documents/GitHub/BrainJ-Python/') # update based on your local file structure

#import QLEAN modules
from BrainJ.Environment import main, imgan

#Settings for Jupyter Notebook
#%matplotlib inline
#%config InlineBackend.figure_format = 'retina'
#%gui qt
# other options
#plt.rcParams["figure.figsize"]=20,20 # makes inline figures full width

#Testing unet initalization and GPU access

main.check_gpu()

#%% Provide Brain dir or dirs: - convert to loop

brain_dirs = [#"Y:/Salzman/Brain-Path]]

for brain_dir in brain_dirs:
    
    # Create a logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    # create a file handler
    logname = brain_dir+'BrainJ_Log_' + datetime.now().strftime("%Y-%m-%d_%H-%M-%S") +'.log'
    handler = logging.FileHandler(logname)
    handler.setLevel(logging.DEBUG)
    
    # create a logging format
    #formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    formatter = logging.Formatter('%(message)s')
    handler.setFormatter(formatter)
    
        # create a stream handler for stdout
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(formatter)
    
    # add the handlers to the logger
    logger.addHandler(handler)
    logger.addHandler(ch)
    #logname = brain_dir+'BrainJ_Log_' + datetime.now().strftime("%Y-%m-%d_%H-%M-%S") +'.log'
    #logging.basicConfig(filename=logname, 
    #                    filemode='w', level=logging.DEBUG, format='%(name)s - %(levelname)s - %(message)s')
    #original_stdout = sys.stdout 
    #sys.stdout = open(logname, 'w')
    
    logger.info("BrainJ Version: "+__version__)
    logger.info("Release Date: "+__date__+"\n")    
    logger.info("Processing brain: "+brain_dir+"\n")

    #Load in experiment parameters and analysis settings   
    settings, locations = main.initialize_brainJ(brain_dir)

    #Modify specific parameters and settings:    
    settings.save_intermediate_data = True
    locations.annotations_table = "C:/Users/Alejandro/Desktop/BrainJ/ABA_CCF_25_2017/Atlas_Regions.csv" #location of atlas regions csv
    locations.annotations_image = "C:/Users/Alejandro/Desktop/BrainJ/ABA_CCF_25_2017/Annotation.tif" #location of atlas regions csv

    settings.tiles_for_prediction = (3,2)
    settings.c1_cell_analysis = True
    settings.c2_cell_analysis = True
    settings.c3_cell_analysis = True
    settings.c4_cell_analysis = True
    
    #locations.cell_analysis_out_dir

    #Print out relevant settings:
    logger.info(f"Analyze channels (C1, C2, C3, C4) set to: {settings.c1_cell_analysis}, {settings.c2_cell_analysis}, {settings.c3_cell_analysis}, {settings.c4_cell_analysis}")
    logger.info(f"Measure channels (C1, C2, C3, C4) set to: {settings.c1_measure}, {settings.c2_measure}, {settings.c3_measure}, {settings.c4_measure}\n")
    

    #1: restore and detect cells
    rawcells = imgan.cell_detection(settings, locations, logger)

    #2a: transform cells - use for normal processing - 2b or 2c are for reprocessing.
    transformedcells = imgan.transform_cells_V2(rawcells, settings, locations, logger)
    
    #2b: if there was an error with elastix import raw cells using the function below, then transform them
    #transformedcells = imgan.import_and_transform_raw_cells_V2(settings, locations, logger) 
    
    #2c: or import already created transformed cells for reannotation, for example with a different atlas
    #transformedcells = imgan.import_transformed_cells_V2(settings, locations) 
    
    #2d: or import already created transformed cells from earlier pipeline for reannotation, for example with a different atlas
    #transformedcells = imgan.import_V1_transformed_cells_V2(settings, locations) 
    
    
    #Annotate cells and create tables
    imgan.annotate_all_cells(transformedcells, settings, locations, logger)
    
    #sys.stdout = original_stdout

    for handler in logger.handlers[:]:
            handler.close()
            logger.removeHandler(handler)
    logger.info("--------------------------------------------------------------------")
