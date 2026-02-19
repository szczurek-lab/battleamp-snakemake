#!/bin/bash

# Function to display script usage
usage() {
    echo "Usage: $0 --output-path <output_path> --end <END> [--delete-json] [--convert-csv]"
    echo "  --output-path    Specify the output directory path"
    echo "  --end            Number of peptides to download (from 1 to end)"
    echo "  --delete-json    (Optional) Specify whether to delete all created JSON files"
    exit 1
}

# Default values
delete_json=false

# Parse command line options
while [[ $# -gt 0 ]]; do
    key="$1"
    echo $key
    case $key in
        --output-path)
            output_path=$2
            shift
            shift
            ;;
        --end)
            END=$2
            shift
            shift
            ;;
        --delete-json)
            delete_json=true
            shift
            shift
            ;;
        *)  # Unknown option
            usage
            ;;
    esac
done


# Check if mandatory options are provided
if [ -z "$output_path" ] || [ -z "$END" ]; then
    usage
fi

mkdir -p "$output_path"   # Creates the output directory if it doesn't exist already
mkdir -p "$output_path/json"   # Creates the output directory if it doesn't exist already

config_file="$output_path/urls_config.txt"   # Sets the path for the configuration file


if [ -e $config_file ]
then
    rm $config_file
fi

# Create config
echo "Creating config file"

for ((i=1; i<=END; i++)); do    # Loops from 1 to END
    echo "url = https://dbaasp.org/peptides/$i" >> "$config_file"    # Appends URL and output paths to the config file
    echo "output = $output_path/json/$i.json" >> "$config_file"     # using a specific format
done

echo "Downloading $END peptide entries from DBAASP to $output_path"
# Download peptides one by one
{ # try
  ./src/curl -X 'GET' -H 'accept: application/json' \
  --parallel \
  --parallel-immediate \
  --retry 30 \
  --retry-delay 10 \
  --config "$output_path/urls_config.txt" ||
  curl -X 'GET' -H 'accept: application/json' \
  --parallel \
  --parallel-immediate \
  --retry 30 \
  --retry-delay 10 \
  --config "$output_path/urls_config.txt"
  }

# Delete JSON files if specified
if [ "$delete_json" = true ]; then
    echo "Deleting JSON files"
    echo "$output_path/json/*.json"
    rm -r -f "$output_path"/json/*.json
fi