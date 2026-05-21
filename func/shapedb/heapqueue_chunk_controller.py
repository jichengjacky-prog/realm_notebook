#!/usr/bin/env python3
import argparse
import heapq
import os
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description='Loop over chunks, keep top ligands with heapqueue, and combine results.'
    )
    parser.add_argument('max_ligands', type=int, help='Number of top ligands to keep')
    parser.add_argument('shapedb_data_location', help='Path to the ShapeDB root directory')
    parser.add_argument('min_chunk', type=int, help='Minimum chunk index to process (inclusive)')
    parser.add_argument('max_chunk', type=int, help='Maximum chunk index to process (inclusive)')
    parser.add_argument('--workers', type=int, default=1, help='Number of worker processes to use')
    parser.add_argument('--output-prefix', default='combined', help='Output filename prefix')
    parser.add_argument('--output-dir', default=None, help='Directory for combined output file')
    return parser.parse_args()


def normalize_data_path(path):
    if not path.endswith('/'):
        path = path + '/'
    return path


def chunk_to_path(chunk):
    chunk_str = str(chunk).zfill(5)
    superchunk_str = str(int(chunk_str[:3]))
    return superchunk_str, chunk_str


def process_chunk_range(max_ligands, shapedb_data_dir, start_chunk, end_chunk):
    heap = []
    for chunk in range(start_chunk, end_chunk + 1):
        superchunk_str, chunk_str = chunk_to_path(chunk)
        print(f'Processing chunk {chunk_str} ({chunk})')

        for subchunk in range(10):
            subchunk_str = str(subchunk)
            path = os.path.join(shapedb_data_dir, superchunk_str, chunk_str, f'{chunk_str}_{subchunk_str}_nn_filtered.txt')

            if not os.path.isfile(path):
                print(f'  missing: {path}')
                continue

            try:
                with open(path, 'r') as fh:
                    for line in fh:
                        line = line.strip()
                        if not line:
                            continue

                        fields = line.split()
                        if len(fields) < 2:
                            continue

                        conf_name = fields[0]
                        try:
                            conf_score = float(fields[1]) * -1
                        except ValueError:
                            continue

                        entry = (conf_score, conf_name, chunk_str, subchunk_str)
                        if len(heap) < max_ligands:
                            heap.append(entry)
                            if len(heap) == max_ligands:
                                heapq.heapify(heap)
                        else:
                            if conf_score > heap[0][0]:
                                heapq.heapreplace(heap, entry)
            except Exception as exc:
                print(f'  error reading {path}: {exc}')

    return heap


def process_superchunk(args):
    max_ligands, shapedb_data_dir, superchunk, min_chunk, max_chunk = args
    start_chunk = max(min_chunk, superchunk * 100)
    end_chunk = min(max_chunk, superchunk * 100 + 99)
    if start_chunk > end_chunk:
        return []
    print(f'Worker starting superchunk {superchunk} (chunks {start_chunk}-{end_chunk})')
    return process_chunk_range(max_ligands, shapedb_data_dir, start_chunk, end_chunk)


def combine_heaps(max_ligands, partial_heaps):
    heap = []
    for partial in partial_heaps:
        for entry in partial:
            if len(heap) < max_ligands:
                heap.append(entry)
                if len(heap) == max_ligands:
                    heapq.heapify(heap)
            else:
                if entry[0] > heap[0][0]:
                    heapq.heapreplace(heap, entry)
    return heap


def process_chunks(max_ligands, shapedb_data_dir, min_chunk, max_chunk, workers=1):
    if workers <= 1:
        return process_chunk_range(max_ligands, shapedb_data_dir, min_chunk, max_chunk)

    from multiprocessing import Pool
    superchunk_start = min_chunk // 100
    superchunk_end = max_chunk // 100
    superchunk_args = [
        (max_ligands, shapedb_data_dir, sc, min_chunk, max_chunk)
        for sc in range(superchunk_start, superchunk_end + 1)
    ]

    with Pool(processes=workers) as pool:
        results = pool.map(process_superchunk, superchunk_args)

    return combine_heaps(max_ligands, results)


def write_results(heap, shapedb_data_dir, output_dir, output_prefix, max_ligands, min_chunk, max_chunk):
    output_parent = output_dir or os.path.join(shapedb_data_dir, 'combined_results')
    os.makedirs(output_parent, exist_ok=True)
    output_file = os.path.join(
        output_parent,
        f'{output_prefix}_best_{max_ligands}_chunks_{str(min_chunk).zfill(5)}_{str(max_chunk).zfill(5)}.txt'
    )

    sorted_results = sorted(heap, reverse=True)
    with open(output_file, 'w') as fh:
        for conf_score, conf_name, chunk_str, subchunk_str in sorted_results:
            fh.write(f'{conf_score},{conf_name},{chunk_str},{subchunk_str}\n')

    print(f'Wrote {len(sorted_results)} combined entries to {output_file}')
    return output_file


def main():
    args = parse_args()
    shapedb_data_dir = normalize_data_path(args.shapedb_data_location)
    if not os.path.isdir(shapedb_data_dir):
        print(f'Error: shapedb_data_location does not exist: {shapedb_data_dir}', file=sys.stderr)
        sys.exit(1)

    if args.min_chunk < 0 or args.max_chunk < args.min_chunk:
        print('Error: invalid chunk range', file=sys.stderr)
        sys.exit(1)

    heap = process_chunks(args.max_ligands, shapedb_data_dir, args.min_chunk, args.max_chunk, workers=args.workers)
    write_results(heap, shapedb_data_dir, args.output_dir, args.output_prefix, args.max_ligands, args.min_chunk, args.max_chunk)


if __name__ == '__main__':
    main()
