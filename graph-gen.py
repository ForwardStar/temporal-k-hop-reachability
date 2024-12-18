import os
import time
import threading
import tarfile
import gzip
import sys

headers = {'User-Agent':'Mozilla/5.0 (Windows; U; Windows NT 6.1; en-US; rv:1.9.1.6) Gecko/20091201 Firefox/3.5.6'}

def download(url, path):
    try:
        from pathlib import Path
        from tqdm import tqdm
    except:
        print("Installing dependencies...")
        from pip._internal import main
        main(['install', 'pathlib'])
        main(['install', 'tqdm'])
        from pathlib import Path
        from tqdm import tqdm
    from urllib.request import urlopen, Request
    print("Fetching from", url + "...")
    path = Path(path)
    blocksize = 1024 * 8
    blocknum = 0
    retry_times = 0
    while True:
        try:
            with urlopen(Request(url, headers=headers), timeout=3) as resp:
                total = resp.info().get("content-length", None)
                with tqdm(
                    unit="B",
                    unit_scale=True,
                    miniters=1,
                    unit_divisor=1024,
                    total=total if total is None else int(total),
                ) as t, path.open("wb") as f:
                    block = resp.read(blocksize)
                    while block:
                        f.write(block)
                        blocknum += 1
                        t.update(len(block))
                        block = resp.read(blocksize)
            break
        except KeyboardInterrupt:
            if path.is_file():
                path.unlink()
            raise
        except:
            retry_times += 1
            if retry_times >= 20:
                break
            print("Timed out, retrying...")
    if retry_times >= 20:
        if path.is_file():
            path.unlink()
        raise ConnectionError("bad internet connection, check it and retry.")

def showProcess():
    print(waiting_message, end="  ")
    while is_finished is False:
        print('\b-', end='')
        time.sleep(0.05)
        print('\b\\', end='')
        time.sleep(0.05)
        print('\b|', end='')
        time.sleep(0.05)
        print('\b/', end='')
        time.sleep(0.05)
    if is_finished is True:
        print('\bdone')
    else:
        print('\berror!')
    
def takeThird(triple):
    return triple[2]

def move_data_file(source, graph_type, destination):
    if source.endswith(".txt"):
        source = open(os.path.join('datasets', source), "r")
    else:
        source = open(os.path.join(os.path.join('datasets', source), "out." + source), "r")
    lines = source.readlines()
    output = []
    for line in lines:
        output.append(line)
        if graph_type == 'U':
            line = line.split()
            if len(line) >= 2:
                line[0], line[1] = line[1], line[0]
                output.append(' '.join(line) + "\n")
    destination = open(destination, "w")
    destination.writelines(output)
    destination.close()

def normalize(filename, increasing=True):
    lines = open(filename, "r").readlines()
    contents = list()
    nodes = dict()

    for line in lines:
        # omit the comments
        if '%' in line or '#' in line:
            continue
        
        # omit the multiplicity of edges
        line = line.split()
        if increasing:
            contents.append([line[0], line[1], int(float(line[len(line) - 1]))])
        else:
            contents.append([line[0], line[1], -int(float(line[len(line) - 1]))])
        if line[0] not in nodes:
            nodes[line[0]] = 1
        if line[1] not in nodes:
            nodes[line[1]] = 1
    
    # normalize nodes
    node_id = 0
    for k in nodes.keys():
        nodes[k] = node_id
        node_id += 1

    # normalize timestamps
    contents.sort(key=takeThird)
    contents[0].append(0)
    contents[0][0] = str(nodes[contents[0][0]])
    contents[0][1] = str(nodes[contents[0][1]])
    for i in range(1, len(contents)):
        contents[i][0] = str(nodes[contents[i][0]])
        contents[i][1] = str(nodes[contents[i][1]])
        if contents[i][2] == contents[i - 1][2]:
            contents[i].append(contents[i - 1][3])
        else:
            contents[i].append(contents[i - 1][3] + 1)

    # wrap up
    text = ""
    for line in contents:
        text += line[0] + " " + line[1] + " " + str(line[3]) + "\n"
    open(filename, "w").write(text)

if __name__ == "__main__":
    # download datasets
    DATASETS_URL = [("http://konect.cc/files/download.tsv.contact.tar.bz2", 'U'),
                    ("http://konect.cc/files/download.tsv.mit.tar.bz2", 'U'),
                    ("https://snap.stanford.edu/data/cit-Patents.txt.gz", 'D'),
                    ("http://konect.cc/files/download.tsv.link-dynamic-nlwiki.tar.bz2", 'D'),
                    ("http://konect.cc/files/download.tsv.soc-sign-bitcoinotc.tar.bz2", 'D'),
                    ("https://snap.stanford.edu/data/email-Eu-core-temporal.txt.gz", 'D'),
                    ("https://snap.stanford.edu/data/CollegeMsg.txt.gz", 'D'),
                    ("http://konect.cc/files/download.tsv.digg-friends.tar.bz2", "D"),
                    ("http://konect.cc/files/download.tsv.epinions.tar.bz2", "D")]
    if os.path.isdir("datasets") is False or len(os.listdir("datasets")) < len(DATASETS_URL):
        need_download = False
        if os.path.isdir("datasets") is False:
            os.mkdir("datasets")
        for (url, graph_type) in DATASETS_URL:
            path = os.path.join("datasets", url.split('/')[-1])
            if not os.path.exists(path):
                if (path.split('.')[-1] == "bz2" and not os.path.exists(os.path.join("datasets", path.split('.')[2]))) or \
                    (path.split('.')[-1] == "gz" and not os.path.exists(path.split('.')[0] + '.' + path.split('.')[1])):
                        if not need_download:
                            need_download = True
                            print("Downloading datasets...")
                        download(url, path)

    # extract all datasets
    waiting_message = "Extracting datasets..."
    is_finished = False
    thread_extract_datasets = threading.Thread(target=showProcess)
    thread_extract_datasets.start()
    file_ls = os.listdir("datasets")
    for file in file_ls:
        if file.endswith(".tar.bz2"):
            archive = tarfile.open(os.path.join("datasets", file), "r:bz2")
            archive.extractall("datasets")
            os.remove(os.path.join("datasets", file))
        elif file.endswith(".txt.gz"):
            archive = gzip.GzipFile(os.path.join("datasets", file))
            out = open(os.path.join("datasets", file.split('.')[0] + "." + file.split('.')[1]), "wb")
            out.write(archive.read())
            archive.close()
            os.remove(os.path.join("datasets", file))
    is_finished = True
    thread_extract_datasets.join()

    # select a target graph dataset
    file_ls = os.listdir("datasets")
    count = 1
    print("Datasets:")
    print("0. example")
    for file in file_ls:
        if file.endswith(".txt"):
            file = file.split(".")[0]
        print(str(count) + ".", file)
        count += 1
    user_input = None
    if len(sys.argv) >= 2:
        user_input = sys.argv[1]
    else:
        user_input = input("Select a graph dataset (0-" + str(count - 1) + "): ")

    # move data file
    if user_input.strip() in [str(i) for i in range(count)]:
        waiting_message = 'Copying dataset to "graph.txt"...'
        is_finished = False
        thread_move_data_file = threading.Thread(target=showProcess)
        thread_move_data_file.start()
        if (int(user_input) == 0):
            open("graph.txt", "w").write("5 3 0\n2 6 0\n1 4 1\n4 5 1\n2 4 2\n4 6 2\n3 4 3\n4 7 3\n7 6 3\n 2 7 3")
        else:
            file = file_ls[int(user_input) - 1]
            graph_type = None
            for (url, gt) in DATASETS_URL:
                if file in url:
                    graph_type = gt
                    break
            move_data_file(file, graph_type, "graph.txt")
        is_finished = True
        thread_move_data_file.join()
    else:
        print("Invalid input! Program terminated.")
        exit()
    
    # normalize the graph
    waiting_message = "Normalizing the graph..."
    is_finished = False
    thread_normalize = threading.Thread(target=showProcess)
    thread_normalize.start()
    if int(user_input) != 0:
        if file_ls[int(user_input) - 1] == 'cit-Patents.txt':
            normalize("graph.txt", False)
        else:
            normalize("graph.txt")
    is_finished = True
    thread_normalize.join()