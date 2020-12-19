S1 = "ttACCTt"
S2 = "gATGTct"
S3 = "ttGCCGT"

l = 4
y = length(S1)
Xlmer = String[]
Xp_lmer = String[]
X_Prime = String[]
z = length(Xp_lmer)


#Create all possible lmers from S1 (Ending with array of lmers in Xlmer)

for i in 1: (y-l+1)  #
    lmer = S1[i: i+l-1]
    push!(Xlmer,lmer)
end


##Store remaining strings into an array. Iterate through each string to produce lmers.
Xp_lmer = [S2,S3]

z = lastindex(Xp_lmer[1])
t = length(Xp_lmer)


#find lmers of X_prime  (Ending with array of lmers in X_Prime)
for j in 1:t    
    for i in 1: (z-l+1)
        lmerXp = Xp_lmer[j][i: i+l-1]
        push!(X_Prime,lmerXp)
    end
end

#Setting an array of tuples for pair comparisons C(X,Si)
C = []
#C = Array{Tuple{String,String},1}
#X = Array{Tuple{String,String},1}
#count of X_prime and X lmers 
m = length(X_Prime)
k = length(Xlmer)
z = 1 


#Iterating through all Xlmers (count k) and pairing each Xlmer with all XPrime lmers (count m)
for i in 1:k
    for j in 1:m
        X = (Xlmer[i],X_Prime[j])
        push!(C,X)
       
    end
end
print(C)
println()
print("The first pair is: $(C[1])")
println()
print("Type in each index of C is:  $(typeof(C[1]))")





















