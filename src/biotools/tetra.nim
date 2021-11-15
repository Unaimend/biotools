import bio_seq
import bionim/utils
import std/[tables]
import sequtils
import strutils
import os
import argparse

#write unit test by pidgheon hole principle
#calculate # permuantattions through math, that check if every element is unique
proc repeatedPermutations[T](a: openarray[T], n: int): seq[seq[T]] =
  result = newSeq[seq[T]]()
  if n <= 0: return
  for i in 0 .. a.high:
    if n == 1:
      result.add(@[a[i]])
    else:
      for j in repeatedPermutations(a, n - 1):
        result.add(a[i] & j)

proc perm[T](a: openarray[T], n: int, use: var seq[bool]): seq[seq[T]] =
  result = newSeq[seq[T]]()
  if n <= 0: return
  for i in 0 .. a.high:
    if not use[i]:
      if n == 1:
        result.add(@[a[i]])
      else:
        use[i] = true
        for j in perm(a, n - 1, use):
          result.add(a[i] & j)
        use[i] = false

proc permutations[T](a: openarray[T], n: int): seq[seq[T]] =
  var use = newSeq[bool](a.len)
  perm(a, n, use)



#REVERSE BULLSHIT INCLUDED
proc generateFreqTableRev() : OrderedTable[string, int]=
  ### ttps://forum.nim-lang.org/t/2812 CREDITS
  var test1:seq[seq[char]] = repeatedPermutations("ACTG", 4)
  var test: OrderedTable[string, int]
  for v in test1:
    var c: string = $complement(v.map(toNucleotide)) 
    var v2 = v.join("")
    if test.hasKey(c) or test.hasKey(v2):
      continue
    else:
      test[v2] = 0
  #echo test
  test

proc generateFreqTable() : OrderedTable[string, int]=
  ### ttps://forum.nim-lang.org/t/2812 CREDITS
  var test1:seq[seq[char]] = repeatedPermutations("ACTG", 4)
  var test: OrderedTable[string, int]
  for v in test1:
    test[v.join("")] = 0
  test

proc countTetraNucleotide*(sequence: string):CountTable[string]=
  let tetras = countKmers(sequence, 4)
  tetras

proc countTetraNucleotideReversed*(sequence: string):OrderedTable[string, int]=
  let tetras = countKmers(sequence, 4)
  var freqs = generateFreqTableRev()
  for k,v in tetras:
    if 'N' in k:
      continue
    if k in freqs:
      freqs[k] += v
    else:
      freqs[$(k.toNucleotide.complement)] += v
  freqs 

proc countTetraNucleotideAll*(sequence: string):OrderedTable[string, int]=
  let tetras = countKmers(sequence, 4)
  var freqs = generateFreqTable()
  for k,v in tetras:
    if 'N' in k:
      continue
    freqs[k] = v
  freqs 

proc countAllFreqs*(sequences: seq[string]): Table[string, string]=
  var output: Table[string, string]
  let temp = generateFreqTable()
  for k,v in temp:
    output[k] = $0
  var header_flag = true 
  for seq in sequences:
    var count = countTetraNucleotideAll(seq)
    if header_flag:
      var head_str = ""
      output["id"] = ""
      #head_str.add("id,")
      for k,v in count:
        output["id"].add(k & ",")
      header_flag = false
    output[seq] = ""
    for k,v in count:
      output[seq].add($v & ",")
  output

proc countSliding(whole_file: string, reversed=false): OrderedTable[string, string] =
  let size=3000 
  let amount_slid_windows = int(len(whole_file)/size)
  let amount_left = len(whole_file) mod size 
  #echo amount_left
  #echo amount_slid_windows
  #echo len(whole_file)
  
  var output: OrderedTable[string, string]

  var header_flag=true
  #var test_str = ""
  output["id"] = ""
  for i in countup(0, amount_slid_windows-1, 1):
    let window = whole_file[i*size .. (i+1)*size-1]
    var t: OrderedTable[string, int]
    if reversed:
      t  = countTetraNucleotideReversed(window)
    else:
      t  = countTetraNucleotideAll(window)
    var index = 0
    if header_flag:
      for key, val in t:
        if index == (if reversed: 127 else: 255):
          index = index + 1
          output["id"].add(key)
          index=0
        else:
          index = index + 1
          output["id"].add(key & ",")
      header_flag = false
    output[$i] = ""
    for k,v in t:
      if index ==  (if reversed: 127 else: 255):
        index += 1
        index=0
        output[$i].add($v)
      else:
        index += 1
        output[$i].add($v & ",")

    #test_str.add(window)
  #ADD OPTION OT INCLUDE OVERLAPP
  let overlap_str: string = whole_file[amount_slid_windows*size ..  (amount_slid_windows*size)+amount_left-1]
  output



proc countFreqs(sequences: seq[string]): Table[string, string]=
  var output: Table[string, string]
  var header_flag = true 
  for seq in sequences:
    var count = countTetraNucleotide(seq)
    if header_flag:
      var head_str = ""
      output["id"] = ""
      #head_str.add("id,")
      for k,v in count:
        output["id"].add(k & ",")
      header_flag = false
    output[seq] = ""
    for k,v in count:
      output[seq].add($v & ",")
  output

proc countAllFreqs(sequences: seq[Sequence], normalized:bool=false): OrderedTable[string, string]=
  var output: OrderedTable[string, string]
  var header_flag = true 
  for seq in sequences:
    var count = countTetraNucleotideAll($seq.data)
    if header_flag:
      var head_str = ""
      output["id"] = ""
      var index = 0
      for k,v in count:
        if index == 255:
          index = index + 1
          output["id"].add(k)
        else:
          index = index + 1
          output["id"].add(k & ",")
      header_flag = false
    output[seq.id] = ""
    var tetra_sum = 0
    for v in count.values:
      tetra_sum = tetra_sum +  v
    var index = 0
    for k,v in count:
      # If have the keys in out, so dont rely on order
      if index == 255:
        index = index + 1
        if not normalized:
          output[seq.id].add($v )
        else:
          let normalized_value = v/tetra_sum
          output[seq.id].add($normalized_value)
      else:
        index = index + 1
        if not normalized:
          output[seq.id].add($v & ",")
        else:
          let normalized_value = v/tetra_sum
          output[seq.id].add($normalized_value & ",")
  output

proc countAllFreqsSplit(sequences: seq[Sequence], normalized:bool=false, reversed:bool = false): OrderedTable[tuple[id: string, sub: string], string]=
  var output: OrderedTable[tuple[id: string, sub: string], string] 
  var header_flag = true 
  var sub = 0
  for seq in sequences:
    var count = countSliding($seq.data, reversed)
    var beg = true
    for k,v in count:
      echo v
      if header_flag:
        # Addd header to output
        output[($sub, "contig")] = v[0 .. (len(v)-1) ]
        header_flag = false
        # Skip to next interation so we do not add header twice
        continue
      if beg == true:
        # because v does include a header sequence for every sequence we skip the first dic entry for every sequence which is the header line
        beg = false
        continue
       
      output[($sub, seq.id)] = v
      sub += 1
  output
    



when isMainModule:
  var p = newParser:
    option("-i", "--input", help="Comma separated list of input FASTA files")
    flag("-s", "--sliding", help="Comma separated list of input FASTA files")
    flag("-r", "--reversed", help="Comma separated list of input FASTA files")
    option("-o", "--output", help="Folder where to output")

  let args: seq[string] = commandLineParams()
  var files: seq[string]
  var oFolder: string
  var sliding: bool
  var reversed: bool = false

  try:
    var opts = p.parse(args)
    files = opts.input.split(',')
    sliding = opts.sliding
    oFolder = opts.output
    reversed = opts.reversed
  except ShortCircuit as e:
    if e.flag == "argparse_help":
      echo p.help 
      quit(1)

  echo files
    
  for f in files:
    echo "##########################" & " " & f & "  " & "###################################################"
    let fasta = parseFastaFile(f)
    if not sliding:
      echo "FINISHED PARSING SEQUENCES"
      var freqs =  countAllFreqsSplit(fasta.seqs, false)
      freqs.to_csv(f & ".csv")
    else:
      var whole_file = ""
      for seq in fasta.seqs:
        whole_file.add($(seq.data))
      var output = whole_file.countSliding(reversed)
      if reversed:
        output.to_csv(f & ".reversed.csv")
      else:
        output.to_csv(f & ".csv")
      
      #test_str.add(overlap_str)
      #echo test_str == whole_file






