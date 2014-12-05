# julia FUBAR3D.main.jl <gridInfoFile> <distanceMatrixFile> [threshold]
#ex: julia FUBAR3D.main.jl datasets/rhodopsin.nwk.grid_info datasets/distance.matrix.thresh65.csv
#ex: julia FUBAR3D.main.jl datasets/rhodopsin.nwk.grid_info datasets/distance.matrix.thresh65.csv 20


# find conditionals
gridInfoFile = ARGS[1]
distanceMatrixFile = ARGS[2]
if length(ARGS)>2
  thresh = float64(ARGS[3])
else
  thresh = 10
end

#RUN MCMC
include("isingSim.jl")
 R = readdlm(distanceMatrixFile,' ')

for i in 1:size(R)[1]
  for j in 1:size(R)[2]
    if R[i,j] > thresh
      R[i,j] = 0
    end
  end
end
S = sparse(R)
β = 0.0
isingSim(S,β,gridInfoFile,1e4)



function broken_shit()
  try
        println(1)
  R = readdlm(distanceMatrixFile,' ')
          println(1)

  if size(R)[1]!=size(R)[2]
        println(1)
    error("R is not square")
        println(1)
  end
catch e
        println(1)
  try
            println(1)
    R = readdlm(distanceMatrixFile,',')
        println(1)
    if size(R)[1]!=size(R)[2]
        println(1)
      error("R is not square")
    end
  catch e
            println(1)
    try
      R = readdlm(distanceMatrixFile,'\t')
      if size(R)[1]!=size(R)[2]
        error("R is not square")
      end
    catch e
      println(e)
    end
  end
end
end
