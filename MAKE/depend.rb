#!/bin/env ruby
#
# created by Sebastian Ehlert in March 2018 (v 1.0)
# a script to create dependencies for FORTRAN source code
#
# this script by no means can get all dependencies, you can
# simply exploit all tricks in FORTRAN to make it fail, but
# it will get the most resonable ones.
# Stuff like
# >      include
# >     .'stuff.f'
# or
# >      use
# >     . gbobc
# will make it fail,
# don't even try
# >       i
# >     i  n
# >     n   c
# >     c    l
# >     l     u
# >     u      d
# >     d       e
# >     e'stuff.f'
# which is valid FORTRAN77
#
# put all modules not belonging to your code here:
FINTRINSIC = ['iso_c_binding','iso_fortran_env','blas95','lapack95']

class FDep < Array
   # grep for include statements, also check for comments
   def select_inc
      self.select do |line|
         line = $` if line.match '!'
         line.match /\A\s*[iN][nN][cC][lL][uU][dD][eE]\s/ and \
         not line.match /\A[*cC]/
      end
   end
   # grep for use statements
   def select_mod
      self.select do |line|
         line = $` if line.match '!' if line
         line.match /\A\s*[uU][sS][eE]\s/ and \
         not line.match /\A[*cC!]/
      end
   end
   # extract all filesnames from the selected lines
   def to_file
      self.map do |line| 
         ((line.split)[1].gsub /['",]/, '') 
      end.uniq.select do |line| 
         not FINTRINSIC.include? line.downcase
      end
   end
   # special transformation needed additionally for modules
   def to_objs
      (self.to_file).map { |line| line + '.o' }
   end
end

if ARGV.empty?
   puts 'Usage: ruby depend.rb <file>'
   exit
end

#printf ":"
# get all lines with include or use statements
dep = FDep.new ARGF.readlines
# parse file for modules and include statements
mod = FDep.new dep.select_mod
inc = FDep.new dep.select_inc

unless inc.empty?
   inc.to_file.each do |file|
      printf " $(INCDIR)/%s", file
   end
end
unless mod.empty?
   mod.to_objs.each do |file|
      printf " $(OUTDIR)/%s", file
   end
end
printf "\n"
