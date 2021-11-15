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



proc countSliding(whole_file: string): OrderedTable[string, string] =
  let size=3000
  let amount_slid_windows = int(len(whole_file)/size)

  let amount_left = len(whole_file) mod size 
  var output: OrderedTable[string, string]

  var header_flag=true
  #var test_str = ""
  output["id"] = ""
  for i in countup(0, amount_slid_windows-1, 1):
    let window = whole_file[i*size .. (i+1)*size-1]
    var index = 0
    if header_flag:
      output["id"].add("GC")
      header_flag = false
    var v =  $count_gc_content(window)
    output[$i] =  v
  output


proc countAllFreqsSplit(sequences: seq[Sequence], normalized:bool=false, reversed:bool = false): OrderedTable[tuple[id: string, sub: string], string]=
  var output: OrderedTable[tuple[id: string, sub: string], string] 
  var header_flag = true 
  var sub = 0
  output[($sub, "contig")] = "GC"
  for seq in sequences:
    var count = countSliding($seq.data)
    var beg = true
    for k,v in count:
      if beg == true:
        # because v does include a header sequence for every #
        # sequence we skip the first dic entry for every sequence which is the header line
        beg = false
        continue
      output[($sub, seq.id)] = v
      sub += 1
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
    var whole_file = ""
    
    for seq in fasta.seqs:
      whole_file.add($(seq.data))

    var gc = countAllFreqsSplit(fasta.seqs)
    #(fasta.seqs)
    #var gc = count_sliding(whole_file)
    gc.to_csv(f & ".gc.csv")
