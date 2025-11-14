#!/usr/bin/env python3
"""
Example workflow demonstrating tool composition in ToolUniverse.

This script shows a complete drug discovery workflow:
1. Find disease-associated targets
2. Retrieve protein sequences
3. Get structure predictions
4. Screen compound libraries
5. Calculate drug-likeness properties
"""

from tooluniverse import ToolUniverse


def drug_discovery_workflow(disease_efo_id: str, max_targets: int = 3):
    """
    Execute a drug discovery workflow for a given disease.

    Args:
        disease_efo_id: EFO ID for the disease (e.g., "EFO_0000537" for hypertension)
        max_targets: Maximum number of targets to process
    """
    tu = ToolUniverse()
    tu.load_tools()

    print("=" * 70)
    print("DRUG DISCOVERY WORKFLOW")
    print("=" * 70)

    # Step 1: Find disease-associated targets
    print(f"\nStep 1: Finding targets for disease {disease_efo_id}...")
    targets = tu.run({
        "name": "OpenTargets_get_associated_targets_by_disease_efoId",
        "arguments": {"efoId": disease_efo_id}
    })
    print(f"✓ Found {len(targets)} disease-associated targets")

    # Process top targets
    top_targets = targets[:max_targets]
    print(f"  Processing top {len(top_targets)} targets:")
    for idx, target in enumerate(top_targets, 1):
        print(f"    {idx}. {target.get('target_name', 'Unknown')} ({target.get('uniprot_id', 'N/A')})")

    # Step 2: Get protein sequences
    print(f"\nStep 2: Retrieving protein sequences...")
    sequences = []
    for target in top_targets:
        try:
            seq = tu.run({
                "name": "UniProt_get_sequence",
                "arguments": {"uniprot_id": target['uniprot_id']}
            })
            sequences.append({
                "target": target,
                "sequence": seq
            })
            print(f"  ✓ Retrieved sequence for {target.get('target_name', 'Unknown')}")
        except Exception as e:
            print(f"  ✗ Failed to get sequence: {e}")

    # Step 3: Predict protein structures
    print(f"\nStep 3: Predicting protein structures...")
    structures = []
    for seq_data in sequences:
        try:
            structure = tu.run({
                "name": "AlphaFold_get_structure",
                "arguments": {"uniprot_id": seq_data['target']['uniprot_id']}
            })
            structures.append({
                "target": seq_data['target'],
                "structure": structure
            })
            print(f"  ✓ Predicted structure for {seq_data['target'].get('target_name', 'Unknown')}")
        except Exception as e:
            print(f"  ✗ Failed to predict structure: {e}")

    # Step 4: Find binding sites
    print(f"\nStep 4: Identifying binding sites...")
    binding_sites = []
    for struct_data in structures:
        try:
            sites = tu.run({
                "name": "Fpocket_find_binding_sites",
                "arguments": {"structure": struct_data['structure']}
            })
            binding_sites.append({
                "target": struct_data['target'],
                "sites": sites
            })
            print(f"  ✓ Found {len(sites)} binding sites for {struct_data['target'].get('target_name', 'Unknown')}")
        except Exception as e:
            print(f"  ✗ Failed to find binding sites: {e}")

    # Step 5: Virtual screening (simplified)
    print(f"\nStep 5: Screening compound libraries...")
    all_hits = []
    for site_data in binding_sites:
        for site in site_data['sites'][:1]:  # Top site only
            try:
                compounds = tu.run({
                    "name": "ZINC_virtual_screening",
                    "arguments": {
                        "binding_site": site,
                        "library": "lead-like",
                        "top_n": 10
                    }
                })
                all_hits.extend(compounds)
                print(f"  ✓ Found {len(compounds)} hit compounds for {site_data['target'].get('target_name', 'Unknown')}")
            except Exception as e:
                print(f"  ✗ Screening failed: {e}")

    # Step 6: Calculate drug-likeness
    print(f"\nStep 6: Evaluating drug-likeness...")
    drug_candidates = []
    for compound in all_hits:
        try:
            properties = tu.run({
                "name": "RDKit_calculate_drug_properties",
                "arguments": {"smiles": compound['smiles']}
            })

            if properties.get('lipinski_pass', False):
                drug_candidates.append({
                    "compound": compound,
                    "properties": properties
                })
        except Exception as e:
            print(f"  ✗ Property calculation failed: {e}")

    print(f"\n  ✓ Identified {len(drug_candidates)} drug candidates passing Lipinski's Rule of Five")

    # Summary
    print("\n" + "=" * 70)
    print("WORKFLOW SUMMARY")
    print("=" * 70)
    print(f"Disease targets processed: {len(top_targets)}")
    print(f"Protein structures predicted: {len(structures)}")
    print(f"Binding sites identified: {sum(len(s['sites']) for s in binding_sites)}")
    print(f"Compounds screened: {len(all_hits)}")
    print(f"Drug candidates identified: {len(drug_candidates)}")
    print("=" * 70)

    return drug_candidates


def genomics_workflow(geo_id: str):
    """
    Execute a genomics analysis workflow.

    Args:
        geo_id: GEO dataset ID (e.g., "GSE12345")
    """
    tu = ToolUniverse()
    tu.load_tools()

    print("=" * 70)
    print("GENOMICS ANALYSIS WORKFLOW")
    print("=" * 70)

    # Step 1: Download gene expression data
    print(f"\nStep 1: Downloading dataset {geo_id}...")
    try:
        expression_data = tu.run({
            "name": "GEO_download_dataset",
            "arguments": {"geo_id": geo_id}
        })
        print(f"  ✓ Downloaded expression data")
    except Exception as e:
        print(f"  ✗ Failed: {e}")
        return

    # Step 2: Differential expression analysis
    print(f"\nStep 2: Performing differential expression analysis...")
    try:
        de_genes = tu.run({
            "name": "DESeq2_differential_expression",
            "arguments": {
                "data": expression_data,
                "condition1": "control",
                "condition2": "treated"
            }
        })
        print(f"  ✓ Found {len(de_genes.get('significant_genes', []))} differentially expressed genes")
    except Exception as e:
        print(f"  ✗ Failed: {e}")
        return

    # Step 3: Pathway enrichment
    print(f"\nStep 3: Running pathway enrichment analysis...")
    try:
        pathways = tu.run({
            "name": "KEGG_pathway_enrichment",
            "arguments": {
                "gene_list": de_genes['significant_genes'],
                "organism": "hsa"
            }
        })
        print(f"  ✓ Found {len(pathways)} enriched pathways")
        if pathways:
            print(f"    Top pathway: {pathways[0].get('pathway_name', 'Unknown')}")
    except Exception as e:
        print(f"  ✗ Failed: {e}")

    print("\n" + "=" * 70)


if __name__ == "__main__":
    # Example 1: Drug discovery workflow for hypertension
    print("EXAMPLE 1: Drug Discovery for Hypertension")
    candidates = drug_discovery_workflow("EFO_0000537", max_targets=2)

    print("\n\n")

    # Example 2: Genomics workflow
    print("EXAMPLE 2: Genomics Analysis")
    genomics_workflow("GSE12345")
