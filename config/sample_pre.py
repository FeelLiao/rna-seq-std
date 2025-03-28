from pathlib import Path
import argparse
import sys
import subprocess
import importlib.util


def check_and_install_pkg(package_name):
    # check if the package is installed
    if importlib.util.find_spec(package_name) is None:
        print(f"⚠️ {package_name} was not found, installing...")
        try:
            # use pip to install the package
            subprocess.check_call([
                sys.executable,
                "-m", "pip", "install",
                "--quiet",  # quiet mode
                "--user",  # avoid permission error
                package_name
            ])
            print(f"✅ {package_name} installed successfully")
        except subprocess.CalledProcessError:
            print(f"❌ Install failed, please install yourself:\
                   {sys.executable} -m pip install {package_name}")
            sys.exit(1)


if __name__ == "__main__":
    check_and_install_pkg("pandas")
    import pandas as pd

VERSION = "0.1.1"

# argument parser
parser = argparse.ArgumentParser(
    description="Automatic create sample sheet for this RNA-seq pipeline",
    epilog="version: " + VERSION)
parser.add_argument(
    "-i", "--input", help="Path of your raw data [default: %(default)s]",
    default="../rawdata", required=False)
parser.add_argument(
    "-e", "--extension", help="file extension of raw data \
      [default: %(default)s]",
    default="fastq", required=False)
parser.add_argument(
    "-o", "--output", help="output path of sample sheet \
      [default: %(default)s]",
    default=".", required=False)
args = parser.parse_args()

filepath = Path(args.input)
output = Path(args.output)
fileExtension = args.extension

# check if the input path and output path exist
assert filepath.exists(), "input path does not exist"
assert output.exists(), "output path does not exist"

# detect all given files in the input path
fqFilesIter = filepath.glob("*.{}".format(fileExtension))
fqFiles = [i for i in fqFilesIter]

samples = list(set([str(i).split("/")[-1].split("_")[0] for i in fqFiles]))
end = list(set(
    [str(i).split("/")[-1].split("_")[1].split(".")[0] for i in fqFiles]))

assert len(list(samples)) > 0, "no files found in input path, please check \
  your input path and file extension"

assert len(end) == 2, "Sample reads are not equal to 2, \
  please check your file format."

# create sample sheet with pandas
sample_sheet = pd.DataFrame(samples, columns=["sample"])
sample_sheet["group"] = sample_sheet["sample"].apply(lambda x: x.split("-")[0])
sample_sheet["read1"] = sample_sheet["sample"].apply(
    lambda x: "{}_{}.{}".format(Path(filepath, x), end[0], fileExtension))
sample_sheet["read2"] = sample_sheet["sample"].apply(
    lambda x: "{}_{}.{}".format(Path(filepath, x), end[1], fileExtension))
sample_sheet["extra"] = ""

outputFile = Path(output, "sample_sheet.csv").absolute()

print("📂 Input path: {}".format(filepath.absolute()))
print("📂 Your file format: {}".format(fileExtension))
print("🔍 Found {} samples".format(len(samples)))
print(f"📝 Sample sheet created successfully, save to\n{outputFile}")

# save the sample sheet to the output path
sample_sheet.to_csv(outputFile, index=False)
