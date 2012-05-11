#!/afs/cs.stanford.edu/u/sidbatra/ruby/install/bin/ruby

t = ENV['t']

main_dir = "/afs/cs.stanford.edu/u/sidbatra/mytest"
exe_dir = "#{main_dir}"
script_dir = "#{main_dir}"
exp_dir = "#{main_dir}/exp_files"
out_dir = "#{main_dir}/output"

#puts `#{exe_dir}/main #{main_dir}/hvac_mdp2/hvac.xml 3 #{exp_dir}/acc_#{t}.txt #{exp_dir}/hor_#{t}.txt #{exp_dir}/rew_#{t}.txt > #{out_dir}/#{t}.txt`
puts `#{exe_dir}/main #{main_dir}/hvac_mdp2/hvac.xml 2 #{exp_dir}/acc_#{t}.txt  #{exp_dir}/rew_#{t}.txt > #{out_dir}/#{t}.txt`

