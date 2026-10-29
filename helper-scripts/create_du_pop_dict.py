import csv
import json
from pathlib import Path


def main():
    script_dir = Path(__file__).resolve().parent
    input_path = Path(input("CSV path: ")).expanduser()
    output_path = script_dir / "du_pop_dictionary.json"

    # Expected headers:
    # DU#, Population code, Population name
    # The code accepts the exact header names above and common variants.
    with input_path.open(newline="") as csv_file:
        reader = csv.DictReader(csv_file)

        if reader.fieldnames is None:
            raise ValueError("CSV file has no header row.")

        fieldnames = {name.strip().lower(): name for name in reader.fieldnames}

        du_key = fieldnames.get("du#") or fieldnames.get("du")
        code_key = (
            fieldnames.get("population code")
            or fieldnames.get("code")
            or fieldnames.get("population_code")
        )
        population_key = (
            fieldnames.get("population name")
            or fieldnames.get("population")
            or fieldnames.get("population_name")
        )

        missing = [
            name
            for name, key in [
                ("DU#", du_key),
                ("Code", code_key),
                ("Population", population_key),
            ]
            if key is None
        ]
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(missing)}")

        records = []
        for row in reader:
            value = row[du_key].strip()
            if value:
                records.append(
                    {
                        "DU#": value,
                        "Code": row[code_key].strip(),
                        "Population": row[population_key].strip(),
                    }
                )

    with output_path.open("w") as json_file:
        json.dump(records, json_file, indent=2)
        json_file.write("\n")

    print(f"Saved {len(records)} records to {output_path}")


if __name__ == "__main__":
    main()