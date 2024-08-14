import argparse
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
#import threading
from pylsl import StreamInfo, StreamOutlet
# import matplotlib.cm as cm  # Import the color map

# run this line in terminal:
# python check_trigno4.py 172.16.9.89 1 2 3 4 5 6 7 8 9 10 11 12

try:
    import pytrigno
except ImportError:
    import sys
    sys.path.insert(0, '..')
    import pytrigno


def plot_data_streaming(data_buffer, axes, y_limits, sensor_ids, colors):
    for i, (ax, color) in enumerate(zip(axes, colors)):
        ax.clear()
        ax.plot(data_buffer[i], label=f'Sensor {sensor_ids[i]+1}', color=color)
        ax.legend(loc='upper right')
        if i == len(axes) - 1:
            ax.set_xlabel('Sample')
        else:
            ax.set_xticks([])  # Remove x ticks for all but the bottom subplot
        ax.set_ylabel('mV')
        ax.set_ylim(y_limits)  # Set fixed y-axis limits
    plt.pause(0.01)  # Pause to update the plot

def check_emg(host, sensor_ids, y_limits=(-1, 1), plot_update_freq=1): 
    print(f"Attempting to connect to EMG device at {host}")

    total_channels = max(sensor_ids)
    print({total_channels})
    channel_range = (0, total_channels -1)
    print('Channel range:', channel_range)

    # Convert sensor IDs to 0-based indexing
    sensor_ids = [sensor_id - 1 for sensor_id in sensor_ids]
    print(sensor_ids)

    if hasattr(pytrigno, 'TrignoEMG'):
        # 100 ms of EMG data per read operation (samples_per_read=200 with 2000 Hz sampling rate)
        dev = pytrigno.TrignoEMG(channel_range=channel_range, samples_per_read=400, host=host)

        try:
            # Start the device
            dev.start()

            # Initialize data buffer
            buffer_size = 2000*3  # Size of the sliding window
            data_buffer = np.zeros((len(sensor_ids), buffer_size))

            # Initialize LSL stream info and outlet
            info = StreamInfo('EMG', 'EMG', total_channels, dev.rate, 'float32', 'myuid34234')
            outlet = StreamOutlet(info)

            sensor_ids_to_plot = sensor_ids[:4]
            print(sensor_ids_to_plot)
            # Set up the interactive plot
            plt.ion()
            fig, axes = plt.subplots(len(sensor_ids_to_plot), 1, figsize=(9, 2 * len(sensor_ids_to_plot)))

            if len(sensor_ids_to_plot) == 1:
                axes = [axes]  # Ensure axes is always iterable

            for ax in axes:
                ax.set_ylim(y_limits)  # Set fixed y-axis limits initially

            fig.suptitle('EMG Signals', y=0.95)  # Set the figure title

            # Generate a color map
            # colors = cm.viridis(np.linspace(0, 1, len(sensor_ids)))
            # Predefined set of visible colors
            colors = ['red', 'blue', 'green', 'purple', 'orange', 'brown', 'pink', 'gray', 'olive', 'cyan']
            
            #iteration = 0
            while True:  # Continuous data collection
                data = dev.read()
                # print(data.shape)
                # Filter data to include only the selected sensor IDs
                filtered_data = data[sensor_ids, :]
                
                # Update the data buffer
                data_buffer = np.roll(data_buffer, -filtered_data.shape[1], axis=1)
                data_buffer[:, -filtered_data.shape[1]:] = filtered_data

                # Send data to LSL
                outlet.push_chunk(data.T.tolist())

                # Update the plot every 'plot_update_freq' iterations
                #if iteration % plot_update_freq == 0:
                #plot_data_streaming(data_buffer, axes, y_limits, sensor_ids_to_plot, colors)
                #iteration += 1
                
            dev.stop()
            plt.ioff()
            plt.show()

        except Exception as e:
            print(f"Error during EMG check: {e}")
    else:
        print("pytrigno module has no attribute 'TrignoEMG'")

"""def check_accel(host):
    print(f"Attempting to connect to Accel device at {host}")
    if hasattr(pytrigno, 'TrignoAccel'):
        channel_range = (0, 17)
        dev = pytrigno.TrignoAccel(channel_range=(0, 17), samples_per_read=10, host=host)
        try:
            dev.start()

            # Initialize data buffer
            buffer_size = 148*3  # Size of the sliding window
            data_buffer = np.zeros((channel_range[1] - channel_range[0] + 1, buffer_size))

            # Set up the interactive plot for accel data
            plt.ion()
            fig, ax = plt.subplots(figsize=(10, 6))

            while True:  # Continuous data collection
                data = dev.read()
                print(f"Received data: {data.shape}")

                # Update the data buffer
                data_buffer = np.roll(data_buffer, -data.shape[1], axis=1)
                data_buffer[:, -data.shape[1]:] = data

                # Plot the data
                plot_data_streaming(data_buffer[:3], ax, 'Accel Data')

            dev.stop()
            plt.ioff()
            plt.show()

        except Exception as e:
            print(f"Error during Accel check: {e}")
    else:
        print("pytrigno module has no attribute 'TrignoAccel'")"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('host', type=str, help="IP address of the machine running TCU. Default is localhost.")
    parser.add_argument('sensor_ids', type=int, nargs='+', help='List of sensor IDs to use')
    parser.add_argument('--plot_update_freq', type=int, default=10, help='Frequency of plot updates')
    args = parser.parse_args()

    check_emg(args.host, args.sensor_ids, plot_update_freq=args.plot_update_freq)
    # check_accel(args.host)

    """emg_thread = threading.Thread(target=check_emg(args.host))
    accel_thread = threading.Thread(target=check_accel(args.host))

    emg_thread.start()
    accel_thread.start()

    emg_thread.join()
    accel_thread.join()
"""