# organoid_boarder_classifier
Static Boarder Classifier based on convex hull and trigonometry


usage:
    usage: classify_boundery_cells.py [-h] [--surpress_plot] [--xlabel XLABEL]
                                  [--ylabel YLABEL] [--arealabel AREALABEL]
                                  [--imagenr IMAGENR] [--colname COLNAME]
                                  input output

Classifies the cells of a 2D slice of an Organiod according if they are a
boarder cell or not

positional arguments:
  input                 absolute path to csv file, or filename if its in the
                        same folder as the .py file.

  output                output path of the filename in obsolute or filename,
                        if filename it gets saved to folder of .py file

optional arguments:
  -h, --help            show this help message and exit
  --surpress_plot, -p   If set, does not create Figures of the Cell
                        classification
  --xlabel XLABEL, -x XLABEL
                        Manually set x label of csv. Default is
                        "Location_Center_X"
  --ylabel YLABEL, -y YLABEL
                        Manually set x label of csv. Default is
                        "Location_Center_Y"
  --arealabel AREALABEL, -a AREALABEL
                        Manually set label the area of a cell is stored in the
                        csv. Default is "AreaShape_Area"
  --imagenr IMAGENR, -i IMAGENR
                        Manually set label the unique identifier for the image
                        is stored in the csv, deafult is "ImageNumber"
  --colname COLNAME, -c COLNAME
                        Specify name of column the boarder classification is
                        saved to, default is "boarder_cell"
