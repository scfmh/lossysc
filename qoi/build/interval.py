import argparse
import numpy as np

def read_binary_file(file_path):
    with open(file_path, 'rb') as file:
        data = np.fromfile(file, dtype=np.float32)
    return data

# def calculate_metrics(original_data, decompressed_data):
#     rela_errors = np.abs(original_data - decompressed_data)
#     rmse = np.sqrt(np.mean((original_data - decompressed_data)**2))
#     max_error = np.max(np.abs(original_data - decompressed_data))
#     avg_original = np.mean(original_data)
#     avg_decompressed = np.mean(decompressed_data)
#     avg_distance = np.abs(avg_original - avg_decompressed)
#     lower_bound_90 = np.percentile(rela_errors, 5)
#     upper_bound_90 = np.percentile(rela_errors, 95)
#     lower_bound_95 = np.percentile(rela_errors, 2.5)
#     upper_bound_95 = np.percentile(rela_errors, 97.5)
#     lower_bound_99 = np.percentile(rela_errors, 0.5)
#     upper_bound_99 = np.percentile(rela_errors, 99.5)
#     return {
#         "90% Range of Relative Errors": (lower_bound_90, upper_bound_90),
#         "95% Range of Relative Errors": (lower_bound_95, upper_bound_95),
#         "99% Range of Relative Errors": (lower_bound_99, upper_bound_99),
#         "RMSE": rmse,
#         "Max Error": max_error,
#         "Avg Distance between Original and Decompressed": avg_distance
#     }

def calculate_metrics(original_data, decompressed_data):
    rela_errors = np.abs(decompressed_data - original_data)
    rmse = np.sqrt(np.mean((original_data - decompressed_data)**2))
    max_error = np.max(np.abs(original_data - decompressed_data))
    avg_original = np.mean(original_data)
    avg_decompressed = np.mean(decompressed_data)
    avg_distance = np.abs(avg_original - avg_decompressed)

    # Calculate dataRange
    min_original = np.min(original_data)
    max_original = np.max(original_data)
    data_range = max_original - min_original

    relative_errors = max_error / data_range

    return {
        "relativeErrors": relative_errors,
        "RMSE": rmse,
        "Max Error": max_error,
        "Avg Distance between Original and Decompressed": avg_distance
    }

def main():
    parser = argparse.ArgumentParser(description='Calculate error metrics between original and decompressed data.')
    parser.add_argument('-o', '--original', type=str, help='Path to the original binary file', required=True)
    parser.add_argument('-d', '--decompressed', type=str, help='Path to the decompressed binary file', required=True)
    args = parser.parse_args()

    original_data = read_binary_file(args.original)
    decompressed_data = read_binary_file(args.decompressed)

    metrics = calculate_metrics(original_data, decompressed_data)

    for metric_name, value in metrics.items():
        print(f"{metric_name}: {value}")

if __name__ == "__main__":
    main()
