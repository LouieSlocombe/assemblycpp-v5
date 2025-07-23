import subprocess
import sys
import os

def run_assembly_tool(inputs_path, exe_path, indices_path):
    total_time = 0
    results = []

    # Paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    summary_path = os.path.join(script_dir, "output_summary.txt")

    # Validate files
    for path, label in [(inputs_path, "inputs"), (exe_path, "executable"), (indices_path, "indices")]:
        if not os.path.isfile(path):
            print(f"Error: {label} file '{path}' does not exist.")
            return

    # Read input file paths
    with open(inputs_path, 'r') as f:
        input_files = [line.strip() for line in f if line.strip()]

    # Read expected assembly indices
    with open(indices_path, 'r') as f:
        expected_indices = [line.strip() for line in f if line.strip()]

    if len(input_files) != len(expected_indices):
        print("Error: Number of inputs does not match number of expected assembly indices.")
        return

    for input_file, expected_index in zip(input_files, expected_indices):
        input_file = "molfiles/" + input_file
        output_file = f"{input_file}Out"
        pathway_file = f"{input_file}Pathway"
        try:
            # Run the executable
            subprocess.run([exe_path, input_file], check=True)

            # Check output file existence
            if not os.path.isfile(output_file):
                msg = f"Output file '{output_file}' not found."
                print("Warning:", msg)
                results.append(msg)
                continue

            with open(output_file, 'r') as out_f:
                lines = out_f.readlines()

                # Check assembly index line
                if not lines:
                    msg = f"Output file '{output_file}' is empty."
                    print("Warning:", msg)
                    results.append(msg)
                    continue

                correct_index_line = f"{input_file} has assembly index: {expected_index}"
                if lines[0].strip() == correct_index_line:
                    index_check = f"[PASS] {input_file}: Correct assembly index {expected_index}"
                else:
                    index_check = f"[FAIL] {input_file}: Expected index {expected_index}, found '{lines[0].strip()}'"

                results.append(index_check)
                print(index_check)

                # Parse time to completion
                time_line = next((l for l in lines if l.startswith("time to completion:")), None)
                if time_line:
                    time_value = int(time_line.split(":")[1].strip())
                    total_time += time_value
                else:
                    results.append(f"Warning: No time info in '{output_file}'")

            # Remove the output file
            os.remove(output_file)
            os.remove(pathway_file)

        except subprocess.CalledProcessError as e:
            error_msg = f"Error: Executable failed on '{input_file}': {e}"
            print(error_msg)
            results.append(error_msg)
        except Exception as e:
            error_msg = f"Unexpected error with '{input_file}': {e}"
            print(error_msg)
            results.append(error_msg)

    total_msg = f"\nTotal time to completion: {total_time}"
    print(total_msg)
    results.append(total_msg)

    # Write all results to output file
    with open(summary_path, 'w') as out_f:
        for line in results:
            out_f.write(line + "\n")

    print(f"\nResults written to: {summary_path}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 script.py <path_to_assembly_exe>")
        sys.exit(1)

    inputs_txt_path = "inputs"
    indices_path = "speed_test_data"
    exe_path = sys.argv[1]

    run_assembly_tool(inputs_txt_path, exe_path, indices_path)
