import argparse
import pandas as pd
import matplotlib.pyplot as plt

# Set the window size for moving average as a global variable
WINDOW_SIZE = 1

def moving_average(data):
    """Calculate the moving average of the given list using the global WINDOW_SIZE."""
    return pd.Series(data).rolling(window=WINDOW_SIZE).mean().tolist()

def visualize_beb(input_file, output_file):
    # Read the table
    df = pd.read_csv(input_file, sep='\t')
    
    # Calculate the moving average for Postmean_W
    smoothed_data = moving_average(df['Postmean_W'])
    
    # Calculate dynamic width based on the number of positions
    DEFAULT_POSITIONS = 100
    BASE_WIDTH = 12
    ratio = len(df) / DEFAULT_POSITIONS
    dynamic_width = BASE_WIDTH * ratio
    
    # Set up the plot with the dynamic width
    fig, ax = plt.subplots(figsize=(dynamic_width, 6))
    
    # Plotting the smoothed Postmean_W values as a simple line plot
    ax.plot(df['Position'], smoothed_data, color='blue', linewidth=1)
    
    ax.set_xlabel("Position")
    ax.set_ylabel("Smoothed Postmean_W")
    ax.set_title("Smoothed Postmean_W vs. Position")
    
    # Set x-axis to start at position 0
    ax.set_xlim(left=0)
    
    # Remove the grid for minimalism
    ax.grid(False)
    
    # Save the plot
    plt.tight_layout()
    plt.savefig(output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize the smoothed Postmean_W from the BEB table.")
    parser.add_argument('-i', '--input', required=True, help="Input table file.")
    parser.add_argument('-o', '--output', required=True, help="Output PNG file for the visualization.")
    args = parser.parse_args()

    visualize_beb(args.input, args.output)
