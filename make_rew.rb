#!/usr/bin/ruby

require 'fileutils'

rew=ARGV[0]
out=ARGV[1]

file = File.open(out,'w');

file.puts "<experiment attribute=\"scale\" type=\"double\" min=\"#{rew}\" max=\"#{rew}\""
file.puts "	    steps=\"1\" interp=\"exponential\">"
file.puts "  <mdp>"
file.puts "    <reward file=\"R_en.txt\">"
file.puts "  </mdp>"
file.puts "</experiment>"

file.close()
