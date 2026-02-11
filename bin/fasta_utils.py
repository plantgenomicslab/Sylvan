import gzip


def readFasta(path: str, sep=" ", index=0) -> dict:
    """Read a FASTA file (plain or gzipped) into a dictionary of {id: sequence}."""
    seq = {}
    seq_id = None
    lines = []

    open_func = gzip.open if path.endswith('.gz') else open
    mode = 'rt' if path.endswith('.gz') else 'r'

    with open_func(path, mode) as file:
        while True:
            line = file.readline()
            if not line:
                break

            line = line.strip()
            if line.startswith('>'):
                if lines:
                    seq[seq_id] = "".join(lines)
                seq_id = line[1:]
                try:
                    seq_id = seq_id.split(sep)[index].strip()
                except (IndexError, ValueError):
                    seq_id = seq_id.split(" ")[0]
                lines = []
            else:
                lines.append(line)

    if seq_id and lines:
        seq[seq_id] = "".join(lines)

    return seq
