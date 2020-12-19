function hamming(x,x1)
    i=1
    count = 0
    while i<=length(x)
        if x[i]!=x1[i]
            count+=1
        end
        i+=1
    end
    return count
end

function filteringcond2(x,x1,z,d)
    i = 1
    count1 = 0
    count2 = 0
    while i<=length(x)
        if x[i]==x1[i]
            if z[i]!=x[i] || z[i]!=x1[i]
                count1+=1
            end
        else
            if z[i]!=x[i]  || z[i]!=x1[i]
                count2+=1
            end
        end
        i+=1
    end
    return [count1,count2]
end

function calculate_alpha_beta(x,x1,z,d)
    alpha_range = length(x)-hamming(x,x1)
    beta_range = hamming(x,x1)
    alphabeta = []

    filtering2 = filteringcond2(x,x1,z,d)

    for i in 0:alpha_range
        for j in 0:beta_range
            k = 2*i+j+hamming(x,x1)
            if k<=2*d
                if abs(filtering2[1]-i)+abs(filtering2[2]-j)<=d
                    return true
                end
            end
        end
    end
    return false
end

print(calculate_alpha_beta("aaggttc","aaccttc","aagctcc",2))