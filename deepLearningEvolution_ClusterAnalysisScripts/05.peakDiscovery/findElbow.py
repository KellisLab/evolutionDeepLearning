import sys
import pandas as pd
import numpy as np

def find_elbow(file_path):
    # Load the cutoff-analysis report
    data = pd.read_csv(file_path, delim_whitespace=True, names=["score", "npeaks", "lpeaks", "avelpeak"])
    data = data[pd.to_numeric(data["score"], errors="coerce").notnull()]
    data["score"] = pd.to_numeric(data["score"])
    data["npeaks"] = pd.to_numeric(data["npeaks"])

    # Get the first and last points
    start_point = data.iloc[0]
    end_point = data.iloc[-1]

    # Calculate the slope and intercept of the line
    slope = (end_point["npeaks"] - start_point["npeaks"]) / (end_point["score"] - start_point["score"])
    intercept = start_point["npeaks"] - slope * start_point["score"]

    # Calculate the distance of each point to the line
    data["distance"] = np.abs(
        slope * data["score"] - data["npeaks"] + intercept
    ) / np.sqrt(slope**2 + 1)

    # Find the index of the maximum distance (elbow point)
    elbow_index = data["distance"].idxmax()

    # Get the corresponding cutoff score
    optimal_cutoff = data.loc[elbow_index, "score"]

    return optimal_cutoff

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python find_elbow.py <cutoff_analysis.txt>")
        sys.exit(1)

    input_file = sys.argv[1]
    cutoff = find_elbow(input_file)
    print(cutoff)
