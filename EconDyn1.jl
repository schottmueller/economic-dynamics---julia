using PyPlot

function ex4122()
    h(x) = 4x*(1-x)
    function quad_ds(start::Float64,n::Int64)
        x = Vector{Float64}(n)
        x[1] = start
        for i = 2:n
            x[i] = h(x[i-1])
        end
        x
    end
    fig,(ax1,ax2,ax3,ax4) = subplots(2,2)
    ax1[:hist](quad_ds(0.1,5000),bins = 40,normed=true,alpha=0.5)
    ax2[:hist](quad_ds(0.2,5000),bins = 40,normed=true,alpha=0.5)
    ax3[:hist](quad_ds(0.3,5000),bins = 40,normed=true,alpha=0.5)
    ax4[:hist](quad_ds(0.6,5000),bins = 40,normed=true,alpha=0.5)
end


function ex4123()
    h(r,x) = r*x*(1-x)
    m = length(2.7:0.01:4.0)
    X = Vector{Float64}(m*50)
    Y = Vector{Float64}(m*50)
    for (j,r) in enumerate(2.7:0.01:4.0)
        x = 0.3
        for i in 1:951
            x = h(r,x)
        end
        Y[(j-1)*50+1] = x
        X[(j-1)*50+1] = r
        for i in (j-1)*50+2:j*50
            Y[i] = h(r,Y[i-1])
            X[i] = r
        end
    end
    #plot(X,Y,"bo",alpha=0.5,linewidth = 0.01)
    scatter(X,Y,s=1,alpha = 0.5)
end


type FinMC
    p::Matrix{Float64}
    x::Int
end

function sample(phi::Vector{Float64})
    a = 0.0
    r = rand()
    for (i,prob) in enumerate(phi)
        if a < r <= a+prob
            return i
        end
        a += prob
    end
end


function FinMCsample(MC::FinMC,n::Int64)
    path = zeros(Int64,n)
    for i in 1:n
        j = sample(vec(MC.p[MC.x,1:end]))
        MC.x = j
        path[i] = j
    end
    return path
end

function list45()
    pH = [0.971 0.029 0.0; 0.145 0.778 0.077; 0.0 0.508 0.492]
    h = FinMC(pH,sample([0.3,0.4,0.3]))
    T1 = FinMCsample(h,1000)
end

function ex423()
    pH = [0.971 0.029 0.0; 0.145 0.778 0.077; 0.0 0.508 0.492]
    h = FinMC(pH,3)
    function marginalFinMC(h::FinMC,phi::Vector{Float64},t::Int64,n::Int64)
        Xt = Vector{Int}(n)
        for j = 1:n
            h.x = sample(phi)
            Xt[j] = FinMCsample(h,t)[end]
        end
        return hist(Xt,size(h.p,2))[2]./n
    end
    marginalFinMC(h,vec([0.0,0.0,1.0]),10,100000)
end

function ex424()
    pH = [0.971 0.029 0.0; 0.145 0.778 0.077; 0.0 0.508 0.492]
    h = FinMC(pH,3)
    function marginalFinMC(h::FinMC,phi::Vector{Float64},t::Int64,n::Int64)
        Xt = zeros(Int64,size(h.p,2))
        for j = 1:n
            h.x = sample(phi)
            Xt[FinMCsample(h,t)[end]] += 1
        end
        return Xt./n
    end
    marginalFinMC(h,vec([0.0,0.0,1.0]),10,100000)
end

function ex4212ff()
    function expPi(h::FinMC,t::Int64,pi::Vector{Float64}, phi::Vector{Float64})
        return phi'*h.p^t *pi
    end
    function expPi(h::FinMC,t::Int64,pi::Vector{Float64})# =ones(Int64,length(pi))  
        return h.p^t *pi
    end
    pH = [0.971 0.029 0.0; 0.145 0.778 0.077; 0.0 0.508 0.492]
    h = FinMC(pH,1)
    println(expPi(h,5,vec([1000.0,0.0,-1000.0]),vec([1.0,0.0,0.0])))
    println(expPi(h,5,vec([1000.0,0.0,-1000.0]),vec([0.0,0.0,1.0])))
    println(expPi(h,1000,vec([1000.0,0.0,-1000.0])))
    println(expPi(h,5,vec([1000.0,0.0,-1000.0]),vec([0.2,0.2,0.6])))
end

function path_prob(p::Matrix{Float64},phi::Vector{Float64},X::Vector{Int64})
        prob = phi[X[1]]
        for t in 2:length(X)
            prob = prob * p[X[t-1],X[t]]
        end
        return prob
    end

function ex4215()
    pH = [0.971 0.029 0.0; 0.145 0.778 0.077; 0.0 0.508 0.492]
    path_prob(pH,[0.2, 0.2, 0.6],[1, 2, 1])
end

function ex4216()
    function path_prob_rec(p::Matrix{Float64},phi::Vector{Float64},X::Vector{Int64},prob::Float64=1.0)
        if length(X)>1
            return path_prob_rec(p,phi,X[1:end-1],prob*p[X[end-1],X[end]])
        else
            return prob*phi[X[1]]
        end
    end
    pH = [0.971 0.029 0.0; 0.145 0.778 0.077; 0.0 0.508 0.492]
    path_prob_rec(pH,[0.2, 0.2, 0.6],[1, 2, 1])
end

function  ex4218()
    function set_path_prob(p::Matrix{Float64},phi::Vector{Float64},D::Matrix{Int64})
        function iterD(p::Matrix{Float64},D::Matrix{Int64},x::Int64)
            sum = 0.0
            if size(D)[1]>1
                for i in D[1,1:end]
                    sum+= p[x,i]*iterD(p,D[2:end,1:end],i)
                end
            else
                for i in D[1,1:end]
                    sum+= p[x,i]
                end
            end
            return sum
        end
        sum = 0.0
        for i in D[1,1:end]
            sum += phi[i]*iterD(p,D[2:end,1:end],i)
        end
        return sum
    end
    pH = [0.971 0.029 0.0; 0.145 0.778 0.077; 0.0 0.508 0.492]
    set_path_prob(pH,[0.2,0.2,0.6],[2 3;2 3;2 3])
end

function ex4219()
    pH = [0.971 0.029 0.0; 0.145 0.778 0.077; 0.0 0.508 0.492]
    h = FinMC(pH,sample([0.2,0.2,0.6]))
    count = 0
    for t in 1:10000
        h.x = sample([0.2,0.2,0.6])
        path = FinMCsample(h,3)
        if [path[i]!=1 for i in 1:3]==[true , true ,true]
            count += 1
        end
    end
    return count/10000
end

function ex4220()
    function pi_path(p::Matrix{Float64},path::Vector{Int64},profits::Vector{Float64},rho::Float64)
        sum = 0.0
        for t in 1:length(path)
            sum+=rho^(t-1) * profits[path[t]]
        end
        return sum
    end
    profits = [1000.0,0.0,-1000.0]
    rho = 1.0/(1.05)
    pH = [0.971 0.029 0.0; 0.145 0.778 0.077; 0.0 0.508 0.492]
    paths = [[t1,t2,t3] for t1 in 1:3, t2 in 1:3, t3 in 1:3]
    pi = 0.0
    for path in paths
        pi+= path_prob(pH,[0.2,0.2,0.6],path)*pi_path(pH,path,profits,rho)
    end
    return pi
end

function ex4222()
    profits = [1000.0,0.0,-1000.0]
    rho = 1.0/(1.05)
    pH = [0.971 0.029 0.0; 0.145 0.778 0.077; 0.0 0.508 0.492]
    h = FinMC(pH,sample([0.2,0.2,0.6]))
    pi = 0.0
    for i in 1:10000
        h.x = sample([0.2,0.2,0.6])
        pi += profits[h.x]
        path = FinMCsample(h,2)
        for t in 1:length(path)
            pi += (rho^(t)) * profits[path[t]]
        end
    end
    return pi/10000
end


function ex4223()
    profits = [1000.0,0.0,-1000.0]
    rho = 1.0/(1.05)
    pH = [0.971 0.029 0.0; 0.145 0.778 0.077; 0.0 0.508 0.492]
    function Epi(p::Matrix{Float64},T::Int64,rho::Float64,phi::Vector{Float64},profits::Vector{Float64})
        sum = 0.0
        for t in 1:T
            sum += (rho^(t-1)* phi' * p^(t-1) * profits)[1]
        end
        return sum
    end
    println(Epi(pH,3,1/1.05,[0.2,0.2,0.6],[1000.0,0.0,-1000.0]))
    pit = [Epi(pH,T,1/1.05,[0.2,0.2,0.6],[1000.0,0.0,-1000.0]) for T in 1:100]
    println(findfirst(i->i>=0,pit))
    plot(1:100,pit,alpha=0.5)
end
            
function ex438()
    pQ = [0.97 0.03 0.0 0.0 0.0;0.05 0.92 0.03 0.0 0.0;0.0 0.04 0.92 0.04 0.0;0.0 0.0 0.04 0.94 0.02;0.0 0.0 0.0 0.01 0.99]
    A = (eye(5)-pQ+ones(5,5))'
    b = ones(5)
    x=\(A,b)
    println(x)
    fig, ax = subplots()
    ax[:bar](1:5,x,0.5,color="blue",alpha=0.5)
end

function ex439()
    pH = [0.971 0.029 0.0; 0.145 0.778 0.077; 0.0 0.508 0.492]
    A = (eye(3)-pH+ones(3,3))'
    b = ones(3)
    x=\(A,b)
    println(x)
    profits = x'*[1000,0,-1000]
end

function ex4310()
    pH = [0.971 0.029 0.0; 0.145 0.778 0.077; 0.0 0.508 0.492]
    A = (eye(3)-pH+ones(3,3))'
    b = ones(3)
    phi=\(A,b)
    h = FinMC(pH,sample(phi))
    X = Array{Int64}(10000)
    for i in 1:10000
        h.x = sample(phi)
        X[i] = FinMCsample(h,20)[end]
    end
    println("simulated dist at T=20: ",hist(X,3)[2]/10000," true dist: ",phi)
end

function ex4325()
    phi = [1.0,0.0,0.0,0.0,0.0]'
    pQ = [0.97 0.03 0.0 0.0 0.0;0.05 0.92 0.03 0.0 0.0;0.0 0.04 0.92 0.04 0.0;0.0 0.0 0.04 0.94 0.02;0.0 0.0 0.0 0.01 0.99]
    phi_new = phi * pQ
    while sum([abs(phi_new[i]-phi[i]) for i in 1:length(phi)])>0.0001
        phi = phi_new
        phi_new = phi * pQ
    end
    return phi_new
end

function ex4326()
    function dobrushin(p::Matrix{Float64})
        d = 1.0
        n,m = size(p)
        for i in 1:n
            for j in 1:n
                a = sum([min(p[i,k],p[j,k]) for k in 1:m])
                if a<d
                    d = a
                end
            end
        end
        return d
    end
    function repeat_dobrushin(p::Matrix{Float64},T::Int64=100)
        q = copy(p)
        for t in 1:T
            a = dobrushin(q)
            if a>0
                return a,t
            else
                q = q*p
            end
        end
        println("Dobrushin 0 for all t up to ",T)
    end
    pQ = [0.97 0.03 0.0 0.0 0.0;0.05 0.92 0.03 0.0 0.0;0.0 0.04 0.92 0.04 0.0;0.0 0.0 0.04 0.94 0.02;0.0 0.0 0.0 0.01 0.99]
    repeat_dobrushin(pQ)
end

function ex4327()
    function meanret(h::FinMC)
        #n,m = size(h.p)
        start = copy(h.x)
        X = Vector{Int64}(10000)
        for i in 1:10000
            h.x = sample(vec(h.p[start,1:end]))
            t = 1
            while h.x != start
                h.x = sample(vec(h.p[h.x,1:end]))
                t += 1
            end
            X[i] = t
        end
        return mean(X)
    end
    pH = [0.971 0.029 0.0; 0.145 0.778 0.077; 0.0 0.508 0.492]
    for i in 1:3
        h = FinMC(pH,i)
        println(1/meanret(h))
    end
end

function ex4329()
    function phi(Q=5,q=2)
        pQ = zeros(Q+1,Q+1)
        for i in 1:q+1
            pQ[i,2:Q+1] = [0.5^(Q+2-j) for j in 2:Q+1]
            pQ[i,1] = 1.0-sum(pQ[i,2:end])
        end
        for i in q+2:Q+1
            pQ[i,2:i] = [0.5^(i+1-j) for j in 2:i]
            pQ[i,1] = 1.0-sum(pQ[i,1:end])
        end
        A = (eye(Q+1)-pQ+ones(Q+1,Q+1))'
        b = ones(Q+1)
        phi = \(A,b)
        return phi
    end
    phi()
end

function ex4330()
    function Epi(q,Q=20,C=0.1)
        pQ = zeros(Q+1,Q+1)
        for i in 1:q+1
            pQ[i,2:Q+1] = [0.5^(Q+2-j) for j in 2:Q+1]
            pQ[i,1] = 1.0-sum(pQ[i,2:end])
        end
        for i in q+2:Q+1
            pQ[i,2:i] = [0.5^(i+1-j) for j in 2:i]
            pQ[i,1] = 1.0-sum(pQ[i,1:end])
        end
        A = (eye(Q+1)-pQ+ones(Q+1,Q+1))'
        b = ones(Q+1)
        phi = \(A,b)
        function pi(X,D)
            order = 0.0
            if X<=q
                order = 1.0
            end
            return min(D,X+(Q-X)*order)-C*order
        end
        g(X) = sum( [pi(X,D)/(2^(D+1)) for D in 0:30] )
        exp_pi = sum([phi[i]*g(i-1) for i in 1:Q+1])
        return exp_pi
    end
    pistart = 0.0
    qopt = 0
    for q in 0:20
        profit =  Epi(q)
        if profit>pistart
            pistart = copy(profit)
            qopt = copy(q)
        end
    end
    return qopt, pistart    
end

function ex4334()
    pH = [0.971 0.029 0.0; 0.145 0.778 0.077; 0.0 0.508 0.492]
    h = FinMC(pH,3)
    path1 = FinMCsample(h,50000)
    return hist(path1)[2]/50000
end

function ex4336()
    pH = [0.971 0.029 0.0; 0.145 0.778 0.077; 0.0 0.508 0.492]
    h = FinMC(pH,3)
    path1 = FinMCsample(h,50000)
    profits = [1000,0,-1000]
    return sum([profits[path1[i]] for i in 1:50000])/50000
end

