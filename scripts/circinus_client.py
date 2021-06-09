#! /usr/bin/env python3

import atexit
import readline
import struct
import sys
import os.path as osp
import glob as gb
from argparse import ArgumentParser

import zmq


def parse_args():
  parser = ArgumentParser("Circinus Client")
  parser.add_argument(
    "-g",
    "--graph",
    help=
    "The name of the graph to use. If the graph is not yet loaded, the server will load the graph if it is preinstalled in the data_dir. Otherwise an error will occur."
  )
  parser.add_argument("--graph_path", help="The path of the graph to use.")
  parser.add_argument("--load_config", default="", help="Load Configuration")
  parser.add_argument("-q", "--query", help="The path to the query. -g must also be specified for query.")
  parser.add_argument("-c", "--query_config", help="The configuration of query, example of query config: cps=cfl,mo=1 2 3 4,cs=dynamic,limit=100000")
  parser.add_argument("-t", "--shutdown", action="store_true")
  parser.add_argument("-s", "--server_host", default="localhost", help="The hostname for Circinus Server.")
  parser.add_argument("-p", "--server_port", type=int, default=55954, help="The port at which Circinus Server is serving.")
  return parser.parse_args()


def connect_to(hostname, port):
  context = zmq.Context()
  send_sock = zmq.Socket(context, zmq.PUSH)
  send_url = "tcp://{0}:{1}".format(hostname, port)

  print("Connecting to Circinus Server at {}\n".format(send_url))
  send_sock.connect(send_url)
  return send_sock


def send_query(args, send_sock):
  """   msg: graph name, query path, query config
        example of query config: cps=cfl,mo=1 2 3 4,cs=dynamic,limit=100000
    """
  cmds = [pack(args.graph, mode='str'), pack(args.query, mode='str'), pack(args.query_config, mode='str')]
  send_sock.send_multipart(cmds)


def load_graph(args, send_sock):
  cmds = [pack(args.graph_path, mode='str'), pack(args.graph, mode='str'), pack(args.load_config, mode='str')]
  send_sock.send_multipart(cmds)


def shutdown_server(args, send_sock):
  cmds = [pack("exit", mode='str')]
  send_sock.send_multipart(cmds)


def pack(data, mode='int'):
  if mode == 'int':
    return struct.pack('i', data)
  elif mode == 'str':
    data = str(data).encode('utf-8')
    return struct.pack('{0}s'.format(len(data)), data)
  elif mode == 'dummy':
    return struct.pack('x', '')


class CircinusCommandCompleter:

  def __init__(self, send_sock, histfile=osp.join(osp.expanduser("~"), ".circinus_history"), history_length=1000):
    self.init_history(histfile, history_length)
    self.send_sock = send_sock
    self.usage = dict()
    self.usage = {
      "load": ["load <pre_installed_graph_name>", "load <graph_path> <graph_name> [<load_config>]"],
      "query": ["query <graph_name> <query_path> <query_config_kvs>"],
      "shutdown": ["shutdown"]
    }

  def init_history(self, histfile, history_length):
    readline.parse_and_bind("tab: complete")
    if hasattr(readline, "read_history_file"):
      try:
        readline.read_history_file(histfile)
      except FileNotFoundError:
        pass
      atexit.register(self.save_history, histfile, history_length)
    readline.set_completer(self.completer)
    readline.set_completer_delims(' \t\n=;')

  def save_history(self, histfile, history_length):
    readline.set_history_length(history_length)
    readline.write_history_file(histfile)

  def completer(self, text, state):
    # file name
    if osp.isdir(text):
      return gb.glob(osp.join(text, '*'))[state]
    elif text.startswith('/'):
      return gb.glob(text + '*')[state]

    options = [i for i in self.usage.keys() if i.startswith(text)]
    if state < len(options):
      return options[state]
    else:
      return None

  def exec_command(self):
    last_cmd = ""
    while True:
      cmd_line = input("Circinus> ")
      cmd_line = cmd_line.split(";")
      for idx, cmd in enumerate(cmd_line):
        if idx == len(cmd_line) - 1:
          last_cmd = last_cmd + " " + cmd
          continue

        if idx == 0:
          cmd = last_cmd + " " + cmd

        if cmd == "" or len(cmd_line) == 1:
          continue

        cmds = [pack(x, mode='str') for x in cmd.split(" ") if x != ""]
        self.send_sock.send_multipart(cmds)
        last_cmd = ""


if __name__ == '__main__':
  args = parse_args()
  conn = connect_to(args.server_host, args.server_port)
  if args.query is not None:
    send_query(args, conn)
  elif args.shutdown:
    shutdown_server(args, conn)
  else:  # interactive command line
    completer = CircinusCommandCompleter(conn)
    completer.exec_command()
