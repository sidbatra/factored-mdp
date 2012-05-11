#!/usr/bin/ruby

main_dir = "/afs/cs.stanford.edu/u/sidbatra/mytest"
scripts_dir = "#{main_dir}"
exp_dir = "#{main_dir}/exp_files"


acc_sign=[1,-1,-1,1]
acc_step=(1-0.95)/5.0
t=0

h=3
while h <= 3 do
    s=-1e-7
    while s >= -5e-5 do
        a=1
        acc=[0.95,0.05,0.05,0.95]
        while a <= 6 do
            
            `#{scripts_dir}/make_acc.rb #{acc.join(",")} #{exp_dir}/acc_#{t}.txt`
            `#{scripts_dir}/make_hor.rb #{h} #{exp_dir}/hor_#{t}.txt`
            `#{scripts_dir}/make_rew.rb #{s} #{exp_dir}/rew_#{t}.txt`
            
            `qsub #{scripts_dir}/cluster.rb -q quicksail -v t=#{t}` #-v h=#{heating[u]},c=#{cooling[u]},v=#{hvac[u]},e=#{equip[u]},i=#{u}

            acc.each_index{|i| acc[i]=acc[i]+acc_sign[i]*acc_step}     
            a+=1
            t+=1
        end
       s*=1.8 
    end
    h+=1

    exit
end

puts t





