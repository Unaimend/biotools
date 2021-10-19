import bio_seq
import bionim/utils
import std/[tables]
import sequtils
import strutils
import math
import os
import argparse

proc count_gc_content*(sequence: string):float=
  var count: int = 0
  for c in sequence:
    if c == 'G' or c == 'C':
      count = count + 1
  return count/len(sequence)


proc gc_content(sequences: seq[Sequence]): OrderedTable[string, string]=
  var output: OrderedTable[string, string]
  var header_flag = true 
  for seq in sequences:
    var count = count_gc_content($seq.data)
    if header_flag:
      var head_str = ""
      output["id"] = ""
      output["id"].add("GC"  )
      header_flag = false
    output[seq.id] = ""
    output[seq.id].add($count)

    #str.add($len(seq.data))
    #f.writeline(str)
  #f.close() 
  output

when isMainModule:
  var p = newParser:
    option("-i", "--input", help="Comma separated list of input FASTA files")
    option("-o", "--output", help="Folder where to output")

  let args: seq[string] = commandLineParams()
  var files: seq[string]
  var oFolder: string

  try:
    var opts = p.parse(args)
    files = opts.input.split(',')
    oFolder = opts.output
  except ShortCircuit as e:
    if e.flag == "argparse_help":
      echo p.help 
      quit(1)

  echo files

  for f in files:
    echo "##########################" & " " & f & "  " & "###################################################"
    let fasta = parseFastaFile(f)
    echo "FINISHED PARSING SEQUENCES"
    var gc = gc_content(fasta.seqs)
    gc.to_csv(f & "gc.csv")
