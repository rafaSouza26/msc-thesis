#!/usr/bin/env python3

import os
import sys
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

def rds_to_csv(input_file, output_file=None):
    """
    Convert an R .rds file to a CSV file
    
    Parameters:
    input_file (str): Path to the input .rds file
    output_file (str, optional): Path to the output .csv file. If None, will use the same name as input but with .csv extension
    
    Returns:
    str: Path to the created CSV file
    """
    # Check if input file exists
    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
    
    # Generate output filename if not provided
    if output_file is None:
        output_file = os.path.splitext(input_file)[0] + '.csv'
    
    # Set up R environment and load base package
    base = importr('base')
    
    # Activate pandas conversion for R objects
    pandas2ri.activate()
    
    try:
        # Read the RDS file using R's readRDS function
        r_data = base.readRDS(input_file)
        
        # Convert R object to pandas DataFrame
        if isinstance(r_data, robjects.vectors.DataFrame):
            pd_df = pandas2ri.rpy2py(r_data)
        else:
            # Try to handle other R objects
            try:
                pd_df = pandas2ri.rpy2py(r_data)
            except Exception as e:
                raise TypeError(f"Could not convert R object to pandas DataFrame: {e}")
        
        # Write to CSV
        pd_df.to_csv(output_file, index=False)
        print(f"Successfully converted {input_file} to {output_file}")
        return output_file
    
    except Exception as e:
        print(f"Error converting {input_file}: {e}")
        return None

def main():
    # Check command line arguments
    if len(sys.argv) < 2:
        print("Usage: python rds_to_csv.py input.rds [output.csv]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    try:
        rds_to_csv(input_file, output_file)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()