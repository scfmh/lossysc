import numpy as np
def generate_newnode(sample_mean, sample_std, sample_hitratio, full_number):
    generated_values = np.random.normal(loc=sample_mean, scale=sample_std, size=int(full_number* sample_hitratio)).astype(int)
    # print("Generated values:", generated_values)
    unique_values, counts = np.unique(generated_values, return_counts=True)
    node_full = len(unique_values)
    # print("Unique values:", unique_values)
    # print("Counts:", counts, node_full)
    return node_full

# generate_newnode(512, 3.2, 0.98, 10000)