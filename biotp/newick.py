import re
from Bio import Phylo
from collections import defaultdict


def parse_mso_table(mso_file: str) -> dict:
    """
    Parse MSO table (tab-delimited).

    Args:
        mso_file: MSO file.

    Returns:
        dict: {MSO_group: set of gene IDs}
    """
    mso_dict = defaultdict(set)
    with open(mso_file, mode="r") as mso_handle:
        for line in mso_handle:
            li = line.strip().split()
            if li[0].startswith("#") or not li:
                continue
            mso_name = li[0]
            gene_ids = li[1:]
            mso_dict[mso_name].update(gene_ids)
    print(f"[INFO] Processed MSO groups: {len(mso_dict)} with total genes: {sum(len(genes) for genes in mso_dict.values())}")
    return mso_dict


def parse_spo_table(spo_file: str) -> dict:
    """
    Parse SPO table (mixed tab and comma delimited).

    Args:
        spo_file: SPO file.

    Returns:
        dict: {SPO_group: set of gene IDs}
    """
    spo_dict = defaultdict(set)
    with open(spo_file, mode="r") as spo_handle:
        for line in spo_handle:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            li = line.split("\t")
            spo_name = f"SP:{li[0]}"
            for l in li[1:]:
                gene_ids = l.split(",")
                gene_ids = {gene.strip() for gene in gene_ids}
                spo_dict[spo_name].update(gene_ids)
        print(f"[INFO] Processed SPO groups: {len(spo_dict)} with total genes: {sum(len(genes) for genes in spo_dict.values())}")
    return spo_dict


def annotate_tree_mso_spo(input_file: str, output_file: str, mso: str, spo: str) -> None:
    """
    Annotate tree leaves with MSO and SPO groups.

    Args:
        input_file:  Input Newick tree file.
        mso_dict:    MSO file.
        spo_dict:    SPO file.
        output_file: Output Newick tree file.
    """
    tree = Phylo.read(input_file, format="newick")
    mso_dict = parse_mso_table(mso)
    spo_dict = parse_spo_table(spo)

    for leaf in tree.get_terminals():
        leaf_name = leaf.name
        leaf_name = leaf_name.replace('"', '')

        ms_hits = [ms for ms, genes in mso_dict.items() if leaf_name in genes]
        sp_hits = [sp for sp, genes in spo_dict.items() if leaf_name in genes]

        mso_label = ms_hits[0] if ms_hits else "MS:NA"
        spo_label = sp_hits[0] if sp_hits else "SP:NA"

        leaf.name = f"{leaf_name} | {mso_label} | {spo_label}"

    Phylo.write(tree, output_file, format="newick")
    print(f"[INFO] Annotated tree saved to {output_file}")


def clean_leaf_names(tree_str: str) -> str:
    """
    Clean leaf names in a Newick tree string by replacing spaces with underscores.
    Args:
        tree_str: Newick tree string.
    """
    pattern = r'([\w\s]+?)\s*-\s*\d+:\d+\[&&NHX:[^\]]+\]'

    def replacer(match):
        species = match.group(1)
        return species.replace(" ", "_")

    return re.sub(pattern, replacer, tree_str)
