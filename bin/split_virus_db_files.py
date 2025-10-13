#!/usr/bin/env python3

# Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
# This file is part of ViMOP and is licensed under the MIT License.
# See the LICENSE file in the root of this repository for full license details.

import sys
import argparse
from pathlib import Path
from collections import defaultdict

import yaml
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

HEADER_PARTS = 6
ALLOWED_ORIENTATIONS = {'forward', 'reverse', 'Unknown'}


def die(msg: str, code: int = 2):
    print(f'ERROR: {msg}', file=sys.stderr)
    sys.exit(code)


def parse_args():
    p = argparse.ArgumentParser(
        description='Validate YAML + FASTA (streaming) and route sequences into curated/filter FASTAs with an output YAML manifest.'
    )
    p.add_argument('-y', '--yaml', required=True, help='Input YAML file with curated and filters sections.')
    p.add_argument('-f', '--fasta', required=True, help='Input FASTA file with required header schema.')
    p.add_argument('-o', '--outdir', default='.', help='Output directory (default: current dir).')
    p.add_argument('--version', required=True, help='Version string to write into the output YAML.')
    p.add_argument('--description', required=True, help='Description string to write into the output YAML.')
    p.add_argument('--blast-db', default='blast_db', help='Value for all.blast_db in output YAML (default: blast_db).')
    p.add_argument('--blast-prefix', default='ALL', help='Value for all.blast_prefix in output YAML (default: ALL).')
    p.add_argument('--config-out', default='virus.yaml', help='Output YAML manifest filename (default: manifest.yaml).')
    p.add_argument('--strict', action='store_true',
                   help='Fail on warnings (kept for compatibility; most curated conflicts are caught during YAML validation).')
    return p.parse_args()


def load_yaml(path: Path) -> dict:
    try:
        with open(path, 'r', encoding='utf-8') as fh:
            data = yaml.safe_load(fh) or {}
    except yaml.YAMLError as e:
        die(f'YAML parsing failed for {path}: {e}')
    except FileNotFoundError:
        die(f'YAML file not found: {path}')
    return data


def validate_yaml_structure_and_index(cfg: dict) -> tuple[dict, dict]:
    """
    Returns:
      validated_cfg: {'curated': Dict, 'filters': Dict}
      org_index: dict species -> {'curated': code_or_None, 'filters': set_of_codes}
    Also ensures that no organism appears in more than one curated group.
    """
    if not isinstance(cfg, dict):
        die('Top-level YAML must be a mapping/dict.')

    curated_in = cfg.get('curated', {})
    filters_in = cfg.get('filters', {})

    if not isinstance(curated_in, dict):
        die('YAML "curated" must be a mapping of codes to entries.')
    if not isinstance(filters_in, dict):
        die('YAML "filters" must be a mapping of codes to entries.')

    curated = {}
    filters = {}

    # Build validated curated and track curated membership to detect duplicates
    curated_membership = {}  # species -> curated_code
    for code, entry in curated_in.items():
        if not isinstance(code, str) or not code:
            die('Each curated key must be a non-empty string.')
        if not isinstance(entry, dict):
            die(f'Curated entry "{code}" must be a mapping.')
        name = entry.get('name')
        orgs = entry.get('organisms')
        if not isinstance(name, str) or not name.strip():
            die(f'Curated "{code}" must have a non-empty "name" string.')
        if not isinstance(orgs, list) or not all(isinstance(x, str) and x.strip() for x in orgs) or not orgs:
            die(f'Curated "{code}" must have a non-empty list "organisms" of strings.')
        curated[code] = {'name': name, 'organisms': set(orgs)}
        for sp in orgs:
            if sp in curated_membership:
                prev = curated_membership[sp]
                die(f'Organism "{sp}" appears in multiple curated groups: "{prev}" and "{code}". '
                    f'Each organism may belong to at most one curated group.')
            curated_membership[sp] = code

    # Validate filters
    for code, entry in filters_in.items():
        if not isinstance(code, str) or not code:
            die('Each filter key must be a non-empty string.')
        if not isinstance(entry, dict):
            die(f'Filter entry "{code}" must be a mapping.')
        orgs = entry.get('organisms')
        if not isinstance(orgs, list) or not all(isinstance(x, str) and x.strip() for x in orgs) or not orgs:
            die(f'Filter "{code}" must have a non-empty list "organisms" of strings.')
        filters[code] = {'organisms': set(orgs)}

    # Build organism index for O(1) routing
    org_index = defaultdict(lambda: {'curated': None, 'filters': set()})
    for sp, cur_code in curated_membership.items():
        org_index[sp]['curated'] = cur_code
    for filt_code, entry in filters.items():
        for sp in entry['organisms']:
            org_index[sp]['filters'].add(filt_code)

    return {'curated': curated, 'filters': filters}, org_index


def parse_header_fields(rec: SeqRecord) -> dict:
    """
    Expected header form (example):
    >PQ541162.1 |MAG: Mammarenavirus lassaense isolate ...|Arenaviridae|Mammarenavirus lassaense|reverse|L
    """
    desc = rec.description
    parts = [p.strip() for p in desc.split('|')]
    if len(parts) != HEADER_PARTS:
        die(f'FASTA header for seq "{rec.id}" does not contain exactly {HEADER_PARTS} pipe-separated fields: "{desc}"')

    seq_id = rec.id
    description = parts[1]
    family = parts[2]
    species = parts[3]
    orientation = parts[4]
    segment = parts[5]

    if orientation not in ALLOWED_ORIENTATIONS:
        die(f'Invalid orientation "{orientation}" for seq "{seq_id}". Allowed: {", ".join(sorted(ALLOWED_ORIENTATIONS))}.')

    if not segment or not isinstance(segment, str):
        die(f'Missing/invalid segment for seq "{seq_id}".')

    return {
        'seq_id': seq_id,
        'description': description,
        'family': family,
        'species': species,
        'orientation': orientation,
        'segment': segment,
    }


def open_output_files(outdir: Path, curated_codes, filter_codes):
    outdir.mkdir(parents=True, exist_ok=True)
    curated_handles = {}
    filter_handles = {}
    for code in curated_codes:
        curated_handles[code] = open(outdir / f'{code}.fasta', 'w', encoding='utf-8')
    for code in filter_codes:
        filter_handles[code] = open(outdir / f'{code}.fasta', 'w', encoding='utf-8')
    return curated_handles, filter_handles


def close_output_files(handles: dict):
    for fh in handles.values():
        try:
            fh.close()
        except Exception:
            pass


def write_record(handle, rec: SeqRecord):
    SeqIO.write([rec], handle, 'fasta')


def stream_route(fasta_path: Path, cfg: dict, org_index: dict, outdir: Path):
    curated_codes = list(cfg['curated'].keys())
    filter_codes = list(cfg['filters'].keys())
    curated_handles, filter_handles = open_output_files(outdir, curated_codes, filter_codes)

    curated_segments_seen = {code: set() for code in curated_codes}

    try:
        for rec in SeqIO.parse(str(fasta_path), 'fasta'):
            meta = parse_header_fields(rec)
            sp = meta['species']
            idx_entry = org_index.get(sp)

            if not idx_entry:
                # Species not requested in curated or filters: skip silently
                continue

            # Curated: at most one (enforced during YAML validation)
            cur_code = idx_entry['curated']
            if cur_code is not None:
                if meta['segment'] != 'Unknown':
                    write_record(curated_handles[cur_code], rec)
                    curated_segments_seen[cur_code].add(meta['segment'])
                else:
                    die(f'"{meta["seq_id"]}" is in curated "{cur_code}" but segment is Unknown.')

            # Filters: zero to many
            for filt_code in idx_entry['filters']:
                write_record(filter_handles[filt_code], rec)

    except FileNotFoundError:
        close_output_files(curated_handles)
        close_output_files(filter_handles)
        die(f'FASTA file not found: {fasta_path}')
    except Exception as e:
        close_output_files(curated_handles)
        close_output_files(filter_handles)
        die(f'Failed during FASTA streaming: {e}')

    close_output_files(curated_handles)
    close_output_files(filter_handles)

    curated_paths = {code: f'{code}.fasta' for code in curated_codes}
    filter_paths = {code: f'{code}.fasta' for code in filter_codes}

    return {
        'curated_paths': curated_paths,
        'filter_paths': filter_paths,
        'curated_segments_seen': curated_segments_seen,
    }


def build_manifest(cfg: dict, routing_info: dict, blast_db: str, blast_prefix: str,
                   version: str, description: str, input_fasta_name: str) -> dict:
    manifest = {
        'all': {
            'blast_db': blast_db,
            'blast_prefix': blast_prefix,
            'fasta': input_fasta_name,  # points to input FASTA
        }
    }

    curated_out = {}
    for code, entry in cfg['curated'].items():
        curated_out[code] = {
            'fasta': routing_info['curated_paths'][code],
            'name': entry['name'],
            'organisms': sorted(entry['organisms']),
            'segments': (sorted(routing_info['curated_segments_seen'][code]) if routing_info['curated_segments_seen'][code] else []),
        }
    manifest['curated'] = curated_out

    filters_out = {}
    for code in cfg['filters'].keys():
        filters_out[code] = routing_info['filter_paths'][code]
    manifest['filters'] = filters_out

    manifest['version'] = version
    manifest['description'] = description
    return manifest


def write_manifest(manifest: dict, path: Path):
    with open(path, 'w', encoding='utf-8') as fh:
        yaml.safe_dump(manifest, fh, sort_keys=False, allow_unicode=True)


def main():
    args = parse_args()

    yaml_path = Path(args.yaml)
    fasta_path = Path(args.fasta)
    outdir = Path(args.outdir)

    cfg_raw = load_yaml(yaml_path)
    cfg, org_index = validate_yaml_structure_and_index(cfg_raw)

    routing = stream_route(fasta_path=fasta_path, cfg=cfg, org_index=org_index, outdir=outdir)

    manifest = build_manifest(
        cfg=cfg,
        routing_info=routing,
        blast_db=args.blast_db,
        blast_prefix=args.blast_prefix,
        version=args.version,
        description=args.description,
        input_fasta_name=fasta_path.name
    )
    manifest_path = outdir / args.config_out
    write_manifest(manifest, manifest_path)

    print(f'Wrote curated FASTAs: {", ".join(f"{k}={v}" for k, v in routing["curated_paths"].items())}', file=sys.stderr)
    print(f'Wrote filter FASTAs: {", ".join(f"{k}={v}" for k, v in routing["filter_paths"].items())}', file=sys.stderr)
    print(f'Wrote manifest: {manifest_path}', file=sys.stderr)


if __name__ == '__main__':
    main()
