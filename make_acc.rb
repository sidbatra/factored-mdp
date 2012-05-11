#!/usr/bin/ruby

require 'fileutils'

acc=ARGV[0].split(",").join(" ")
out=ARGV[1]

file = File.open(out,'w');

file.puts "<experiment attribute=\"values\" type=\"vector\" min=\"#{acc} \""
file.puts "	    max=\"#{acc} \" steps=\"#{1}\" interp=\"linear\">"
file.puts "  <mdp>"
file.puts "    <prediction node=\"Occ\">"
file.puts "      <prediction_cpd horizon=\"all\"/>"
file.puts "    </prediction>"
file.puts "  </mdp>"
file.puts "</experiment>"

file.close()
