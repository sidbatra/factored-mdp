#!/usr/bin/ruby

require 'fileutils'

hor=ARGV[0]
out=ARGV[1]

file = File.open(out,'w');

file.puts "<experiment attribute=\"horizon\" type=\"integer\" min=\"#{hor}\" max=\"#{hor}\" steps=\"1\">"
file.puts " <mdp>"
file.puts "      <prediction node=\"Occ\"/>"
file.puts "  </mdp>"
file.puts "</experiment>"

file.close()
