"""Package verified research artifacts; remove private runtime paths in copies.

Original local files are read-only inputs. Every public file has both its
original and distributed SHA-256 recorded. No numerical prediction is edited.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import io
import json
import re
import tarfile

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / 'paper'
PRIVATE_PROJECTS = set()


def digest(data):
    return hashlib.sha256(data).hexdigest()


def clean(text):
    text = text.replace(str(ROOT), '${PROJECT_ROOT}')
    for project in PRIVATE_PROJECTS:
        text = text.replace(project, '${GOOGLE_CLOUD_PROJECT}')
    text = re.sub(r'/work/[^/\s"\x27]+/[^/\s"\x27]+', '${AMR_WORK_ROOT}', text)
    text = re.sub(r'/home/[^/\s"\x27]+', '${USER_HOME}', text)
    if '_umass_edu' in text:
        text = re.sub(r'\b[A-Za-z0-9_]+_umass_edu', '${PRIVATE_ACCOUNT}', text)
    text = re.sub(r'[A-Za-z]:\\+(?:Users|users)\\+[^\\/\s"\x27]+', '${LOCAL_USER_HOME}', text)
    return text


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--documents-only', action='store_true')
    args = parser.parse_args()
    destination = PAPER / 'release'
    destination.mkdir(exist_ok=True)
    if not args.documents_only:
        artifacts = PAPER / 'artifacts'
        records, blobs = [], []
        source_files = []
        for folder in ('data/source/ncbi_ast_snapshot', 'data/source/ncbi_isolate_enrichment/20260821T012604Z-v1', 'data/source/jarbs/20260825T021722Z-v1', 'data/interim'):
            source_files.extend(p for p in (ROOT / folder).glob('*') if p.is_file() and p.suffix in {'.json', '.csv', '.tsv', '.sql'} and 'eucast' not in p.name.lower())
        def collect_projects(value):
            if isinstance(value, dict):
                for key, item in value.items():
                    if key in {'query_project', 'project_id', 'projectId'} and isinstance(item, str) and not item.startswith('ncbi-'):
                        PRIVATE_PROJECTS.add(item)
                    collect_projects(item)
            elif isinstance(value, list):
                for item in value:
                    collect_projects(item)
        for path in source_files:
            if path.suffix == '.json':
                collect_projects(json.loads(path.read_text()))
        for path in sorted(list(artifacts.rglob('*')) + source_files):
            if not path.is_file():
                continue
            relative = str(path.relative_to(artifacts)) if path.is_relative_to(artifacts) else 'source_inputs/' + str(path.relative_to(ROOT))
            # Cleanup transcripts are unrelated to reproducing the scientific runs.
            if 'cleanup' in relative.lower():
                continue
            raw = path.read_bytes()
            binary = path.suffix == '.npz'
            public = raw if binary else clean(raw.decode('utf-8')).encode('utf-8')
            assert b'_umass_edu' not in public and b'/work/' not in public and b'/home/' not in public, relative
            record = {'file': relative, 'original_sha256': digest(raw), 'published_sha256': digest(public),
                      'bytes': len(public), 'redacted': raw != public}
            # Frozen numerical prediction files must remain byte-for-byte unchanged.
            if path.name.startswith(('predictions_', 'scored_predictions_')) and path.suffix == '.csv':
                assert raw == public, ('prediction unexpectedly changed', relative)
            records.append(record)
            if binary:
                blobs.append((relative, path, None))
            else:
                blobs.append((relative, None, public))
        release_manifest = {'format': 1, 'description': 'Private runtime paths redacted in distributed copies. Original hashes retained; numerical predictions unchanged.', 'files': records}
        encoded = (json.dumps(release_manifest, indent=2) + '\n').encode()
        (PAPER / 'provenance/release_manifest.json').write_bytes(encoded)
        blobs.append(('release_manifest.json', None, encoded))
        for label, binary in [('results', False), ('embeddings', True)]:
            target = destination / f'amr-{label}-v1.0.0.tar.gz'
            with tarfile.open(target, 'w:gz', compresslevel=3, dereference=True) as tar:
                for relative, source, data in blobs:
                    if (source is not None) != binary:
                        continue
                    name = 'amr-research-artifacts/' + relative
                    if source is not None:
                        tar.add(source, arcname=name, recursive=False)
                    else:
                        info = tarfile.TarInfo(name)
                        info.size = len(data)
                        info.mode = 0o644
                        tar.addfile(info, io.BytesIO(data))
                for notice in ('THIRD_PARTY_NOTICES.md', 'LICENSE'):
                    tar.add(ROOT / notice, arcname='amr-research-artifacts/' + notice, recursive=False)
            print(target.name, target.stat().st_size, flush=True)
        print('Packaged', len(records), 'files;', sum(r['redacted'] for r in records), 'path-redacted copies', flush=True)

    # Source-only arXiv upload: resolved bibliography, vector figures, supplement.
    with tarfile.open(destination / 'arxiv-source-v1.0.0.tar.gz', 'w:gz') as tar:
        for name in ('manuscript.tex', 'supplement.tex', '00README.json'):
            tar.add(PAPER / name, arcname=name)
        for path in sorted((PAPER / 'figures').glob('*.pdf')):
            tar.add(path, arcname='figures/' + path.name)
    with tarfile.open(destination / 'paper-and-tables-v1.0.0.tar.gz', 'w:gz') as tar:
        for path in sorted(PAPER.rglob('*')):
            if not path.is_file() or set(path.relative_to(PAPER).parts) & {'artifacts', 'release', 'tex-build', '__pycache__'}:
                continue
            tar.add(path, arcname='paper/' + str(path.relative_to(PAPER)))
    lines = []
    for path in sorted(destination.glob('*')):
        if path.is_file() and path.name != 'SHA256SUMS':
            with path.open('rb') as handle:
                lines.append(hashlib.file_digest(handle, 'sha256').hexdigest() + '  ' + path.name)
    (destination / 'SHA256SUMS').write_text('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
