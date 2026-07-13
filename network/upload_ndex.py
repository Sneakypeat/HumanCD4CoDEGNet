#!/usr/bin/env python3
"""Upload HumanCD4CoDEGNet.cx2 to NDEx (www.ndexbio.org).

Your NDEx credentials are read from the environment or prompted at runtime — they are
NEVER written into this file or the repo. Uploads as PRIVATE; flip to PUBLIC in the NDEx
web UI (or set NDEX_VISIBILITY=PUBLIC) once you're ready to share.

    pip install ndex2
    python upload_ndex.py
    # or non-interactively:
    NDEX_USER=you NDEX_PASS=secret python upload_ndex.py
"""
import os, getpass, ndex2
from ndex2.cx2 import RawCX2NetworkFactory

CX2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "HumanCD4CoDEGNet.cx2")
user = os.environ.get("NDEX_USER") or input("NDEx username: ")
pw   = os.environ.get("NDEX_PASS") or getpass.getpass("NDEx password: ")
vis  = os.environ.get("NDEX_VISIBILITY", "PRIVATE").upper()

net = RawCX2NetworkFactory().get_cx2network(CX2)
client = ndex2.client.Ndex2(host="www.ndexbio.org", username=user, password=pw)
uri = client.save_new_cx2_network(net.to_cx2(), visibility=vis)
uuid = uri.rstrip("/").split("/")[-1]
print(f"\nUploaded ({vis}). Network UUID: {uuid}")
print(f"View:  https://www.ndexbio.org/viewer/networks/{uuid}")
if vis != "PUBLIC":
    print("It is PRIVATE for now — open it in NDEx and set visibility to PUBLIC to share.")
