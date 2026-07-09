import gzip


def readFasta(path: str, sep=None, index=0) -> dict:
    """Read a FASTA file (plain or gzipped) into a dictionary of {id: sequence}.

    The record ID is the first whitespace-delimited token of the header. Aligners
    such as miniprot report that same token (e.g. `Target=<id>`), so splitting on
    spaces alone would mis-key headers whose fields are tab-delimited. Pass
    `sep`/`index` to slice the token further, e.g. sep="|", index=1.
    """
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
                header = line[1:].split()
                seq_id = header[0] if header else ""
                if sep is not None:
                    try:
                        seq_id = seq_id.split(sep)[index]
                    except IndexError:
                        pass
                lines = []
            else:
                lines.append(line)

    if seq_id and lines:
        seq[seq_id] = "".join(lines)

    return seq
