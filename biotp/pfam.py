def extract_protein_by_entry(input_file: str, output_file: str, *entry: str) -> None:
    """
    Extract unique protein names corresponding to specific Pfam entries.

    Args:
        input_file:  Input filename.
        output_file: Output filename.
        entry:       Pfam entry to filter.

    """
    MIN_SCORE = 30
    if not entry:
        raise ValueError("[ERROR] No Pfam entry provided.")

    # Initialize a dictionary to hold unique proteins for each entry
    entry2protein = {e: set() for e in entry}

    with open(input_file, mode="r") as input_handle:
        for line in input_handle:
            if line.startswith("#") or not line.strip():
                continue

            columns = line.split()
            if len(columns) < 14:
                continue

            pfam = columns[1]
            protein_name = columns[3]

            try:
                domain_score = float(columns[13])
            except ValueError:
                print(f"[WARNING] Invalid score '{columns[13]}' for protein '{protein_name}'. Skipping.")
                continue

            if pfam in entry2protein and domain_score >= MIN_SCORE:
                entry2protein[pfam].add(protein_name)

    if any(len(proteins) == 0 for proteins in entry2protein.values()):
        print("[WARNING] One or more Pfam entries have zero proteins after filtering.")

    common_proteins = set.intersection(*entry2protein.values())

    with open(output_file, "w") as output_handles:
        output_handles.write("\n".join(sorted(common_proteins)))

    print(f"[INFO] Extracted {len(common_proteins)} unique proteins for Pfam entries: {', '.join(entry)}.")
