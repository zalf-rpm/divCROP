#!/usr/bin/python
# -*- coding: UTF-8

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/. */

# Authors:
# Michael Berg-Mohnicke <michael.berg@zalf.de>
#
# Maintainers:
# Currently maintained by the authors.
#
# This file has been created at the Institute of
# Landscape Systems Analysis at the ZALF.
# Copyright (C: Leibniz Centre for Agricultural Landscape Research (ZALF)

from collections import defaultdict, OrderedDict
import csv
from datetime import datetime
import gc
import json
import numpy as np
import os
from pyproj import CRS, Transformer
import sqlite3
import sys
import timeit
import types
import zmq

import monica_io3
import soil_io3
import monica_run_lib as Mrunlib

PATHS = {
    "mbm-local-remote": {
        "path-to-output-dir": "out/",
        "path-to-csv-output-dir": "csv-out/"
    },
    "remoteConsumer-remoteMonica": {
        "path-to-output-dir": "/out/out/",
        "path-to-csv-output-dir": "/out/csv-out/"
    }
}

def run_consumer(leave_after_finished_run = True, server = {"server": None, "port": None}, shared_id = None):
    "collect data from workers"

    config = {
        "mode": "mbm-local-remote",  
        "port": server["port"] if server["port"] else "7777", 
        "server": server["server"] if server["server"] else "login01.cluster.zalf.de", 
        "timeout": 600000 # 10 minutes
    }

    if len(sys.argv) > 1 and __name__ == "__main__":
        for arg in sys.argv[1:]:
            k,v = arg.split("=")
            if k in config:
                config[k] = v

    paths = PATHS[config["mode"]]

    print("consumer config:", config)

    context = zmq.Context()
    socket = context.socket(zmq.PULL)

    socket.connect("tcp://" + config["server"] + ":" + config["port"])
    socket.RCVTIMEO = config["timeout"]
    leave = False

    #count cols in rows
    def process_message(msg):
        if len(msg["errors"]) > 0:
            print("There were errors in message:", msg, "\nSkipping message!")
            return

        if not hasattr(process_message, "msg_count"):
            process_message.msg_count = defaultdict(lambda: defaultdict(lambda: 0))

        print("received work result ", process_message.received_env_count, " customId: ", str(msg.get("customId", "")))

        custom_id = msg["customId"]
        setup_id = custom_id["setup_id"]
        id = custom_id["id"]
        coord_id = custom_id["coord_id"]
        crop_id_short = custom_id["crop_id_short"]
        
        write_monica_csv(msg, count=len(process_message.msg_count[setup_id].keys()), dir=paths["path-to-output-dir"]+crop_id_short, id=id, coord_id=coord_id)

        process_message.msg_count[setup_id][id] += 1

        return leave

    process_message.received_env_count = 1

    while not leave:
        try:
            msg = socket.recv_json() #encoding="latin-1"
            leave = process_message(msg)
        except zmq.error.Again as _e:
            print('no response from the server (with "timeout"=%d ms) ' % socket.RCVTIMEO)
            return
        except Exception as e:
            print("Exception:", e)
            #continue

    print("exiting run_consumer()")
    #debug_file.close()


def write_monica_csv(result, count, dir, id, coord_id):

    bin_size = 1000
    count_bin = int(count / bin_size)
    dir = dir + "/ids_" + str(count_bin*bin_size) + "-" + str((count_bin+1)*bin_size) 
    if os.path.isdir(dir) and os.path.exists(dir):
        pass
    else:
        try:
            os.makedirs(dir)
        except OSError:
            print("c: Couldn't create dir:", dir, "! Exiting.")
            exit(1)

    with open(dir + "/" + str(coord_id) + ".csv", "w", newline="") as _:
        writer = csv.writer(_, delimiter=",")

        for data_ in result.get("data", []):
            results = data_.get("results", [])
            orig_spec = data_.get("origSpec", "")
            output_ids = data_.get("outputIds", [])

            if len(results) > 0:
                writer.writerow([orig_spec.replace("\"", "")])
                for row in monica_io3.write_output_header_rows(output_ids,
                                                                include_header_row=True,
                                                                include_units_row=False,
                                                                include_time_agg=False):
                    writer.writerow(row)

                for row in monica_io3.write_output(output_ids, results, round_ids={
                    "SowDOY": 1,
                    "HarDOY": 1
                }):
                    writer.writerow(row)

            writer.writerow([])


if __name__ == "__main__":
    run_consumer()


