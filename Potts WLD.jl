using Statistics
function latticesq(ll)
    nn = ll * ll  # 总格点数
    nbor = Array{Int}(undef, 8, nn)  # 第1-4列为最近邻，第5-8列为次近邻
    for s0 = 1:nn
        x0 = mod(s0 - 1, ll)
        y0 = div(s0 - 1, ll)
        
        # 最近邻坐标（右、上、左、下）
        x_right = mod(x0 + 1, ll)
        x_left  = mod(x0 - 1 + ll, ll)
        y_up    = mod(y0 + 1, ll)
        y_down  = mod(y0 - 1 + ll, ll)
        
        # 次近邻坐标（右上、右下、左上、左下）
        x_right_up   = mod(x0 + 1, ll)
        y_right_up   = mod(y0 + 1, ll)
        x_right_down = mod(x0 + 1, ll)
        y_right_down = mod(y0 - 1 + ll, ll)
        x_left_up    = mod(x0 - 1 + ll, ll)
        y_left_up    = mod(y0 + 1, ll)
        x_left_down  = mod(x0 - 1 + ll, ll)
        y_left_down  = mod(y0 - 1 + ll, ll)
        
        # 最近邻索引
        nbor[1, s0] = 1 + x_right + y0 * ll          # 右
        nbor[2, s0] = 1 + x0 + y_up * ll              # 上
        nbor[3, s0] = 1 + x_left + y0 * ll           # 左
        nbor[4, s0] = 1 + x0 + y_down * ll           # 下
        
        # 次近邻索引
        nbor[5, s0] = 1 + x_right_up + y_right_up * ll    # 右上
        nbor[6, s0] = 1 + x_right_down + y_right_down * ll # 右下
        nbor[7, s0] = 1 + x_left_up + y_left_up * ll      # 左上
        nbor[8, s0] = 1 + x_left_down + y_left_down * ll  # 左下
    end
    return nbor
end

 # Constructs the initial random spin configuration

function initspin_jitai(ll::Int)
    nx=ll
    nn=ll*ll
    spin=Array{Int32}(undef,nn)
     for i=1:nx
        for j=1:nx
            ns=j+nx*(i-1)
            spin[ns] = (i % 2 == 1) ? (j % 2 == 1 ?  1 : 3) : (j % 2 == 1 ? 3 : 1)
       #spin[ns]=(i+j)%3+1
       end
    end
return spin
  end

 function initspin(nn::Int)
     spin=Array{Int32}(undef,nn)
     for i=1:nn
        spin[i]=rand(1:1:3)
     end
     return spin
  end




function energy(ll::Int, nbor, spin, J1, J2)
  nx=ll;
  nn=ll*ll;
 current_E = 0.0
for i = 1:nx*nx
                 sum_spin_nbor = 0
                 sum_spin_next_nbor = 0                        
                    for knbor = 1:2                         #最近邻
                             id_nbor = nbor[knbor,i]
                                 if spin[id_nbor] == spin[i]
                                      sum_spin_nbor += 1
                                 end
                    end                
                    for knbor = 5:6                          #次近邻
                             id_nbor = nbor[knbor,i]
                                 if spin[id_nbor] == spin[i]
                                     sum_spin_next_nbor += 1
                                 end
                    end                  
             current_E += -J1 * sum_spin_nbor - J2 * sum_spin_next_nbor    # 能量总和
                end
    return current_E
end
     



function Energies(ll::Int)              # 加次近邻的能量   
nx =ll;nn =ll*ll;
number = (nx - 1) * (nx + 1) ÷ 3      #圈数
ii = 4 * nn +1 - 3 - number - 1 - 3    #能量数
Energy = zeros(Int32, ii)  
E_min = -2 * nn;    E_max = 2 * nn;
Energy[1]  = E_min; Energy[2]   = E_min+4;
Energy[ii] = E_max; Energy[ii-1]= E_max-4; Energy[ii-2]=E_max-4-2;
for i in 1:number
    Energy[i + 2] = E_min+4 + 2 * i 
end 
number_2=ii-4-number      #剩余变化为1的能量数
number_3=2+number         #索引
for j in 1:number_2
    Energy[j+number_3]=E_min+j+4+2*number
end
return Energy
end 
    

function Energies6(ll::Int)
    nx = ll
    nn = nx * nx
 #   number = (3 * nn - 4 * nx + 4) ÷ 8  # 使用整数除法确保结果为整数
    number=5
    ii = 4 * nn - 6 - number  # 计算总能量数
    Energy = zeros(Int32, ii)   
    E_min = -2 * nn
    E_max = 2 * nn 
    # 初始化首尾元素
    Energy[1] = E_min
    Energy[2] = E_min + 4
    Energy[ii] = E_max
    Energy[ii-1] = E_max - 4
    Energy[ii-2] = E_max - 6
    
    # 填充第一部分能量
    for i in 1:number
        Energy[i + 2] = E_min + 4 + 2 * i
    end
    

    number_2 = ii - 5 - number  
    number_3 = 2 + number
    
    # 填充第二部分能量
    for j in 1:number_2
        Energy[j + number_3] = E_min + j + 4 + 2 * number
    end
    
    return Energy
end 
    
function isflat(Hist,step)
    # Initialize variables
    flat_thres = 0.9
    oknum = 0
    H_number = length(Hist) # Assuming his array has length 2*Nb + 1
    div_number=H_number
    for i in 1:H_number
        expected=step / div_number
        prob= Hist[i] / expected
        if prob >=flat_thres
            oknum+=1
        end
    end
  #  print("oknum:",oknum,'\n')
    isflat = 0
    if oknum >= H_number
        isflat = 1
    end
end


function find_energy_index(E::Float64, Energy::Vector{Int32})
    idx = searchsortedfirst(Energy, E)
    if idx == 1
        return 1
    elseif idx > length(Energy)
        return length(Energy)
    else
        return abs(E - Energy[idx-1]) < abs(E - Energy[idx]) ? idx-1 : idx
    end
end 

function Wanglandau(ll::Int,J1,J2,nbor,spin,mstep)
        nx=ll
        nn=ll*ll
        qn = 3  
Energy = Energies6(ll)
 #   print(Energy,'\n')
#Ene=-2*nn
    spin_jitai=initspin_jitai(nx)
Ene=energy(nx, nbor, spin_jitai, J1, J2)
    print(Ene,'\n')
E_number = length(Energy)                  # 定义能量区间               
      S  = zeros(Float64,E_number)       # 态密度G(E)取对数，初始化为0      G(E)=exp^(S(E))  
   Hist  = zeros(Float64,E_number)       # 直方图记录访问次数  
#G(E)=exp^(S(E))  初始化,S(E)=0 G(E)=1
 f=1.0; step=0;
for m=1:mstep
ns=rand(1:nn)
noldspin = spin[ns];
allowed_spins = setdiff(1:3, spin[ns]) 
t_value = rand(allowed_spins)             # 翻转后的自旋 
sum_spin_nbor_old = 0;sum_spin_next_nbor_old = 0;  
sum_spin_nbor_new = 0;sum_spin_next_nbor_new = 0;
   for knbor = 1:4                         #最近邻
        id_nbor = nbor[knbor,ns]
                if spin[id_nbor] == noldspin
                      sum_spin_nbor_old  += 1
                elseif spin[id_nbor] == t_value
                       sum_spin_nbor_new += 1
                end
    end                
    for knbor = 5:8                          #次近邻
        id_nbor = nbor[knbor,ns]
                if spin[id_nbor] == noldspin
                      sum_spin_next_nbor_old += 1
                elseif spin[id_nbor] == t_value
                      sum_spin_next_nbor_new += 1 
                end
    end
    E_total_old=-J1*sum_spin_nbor_old-J2*sum_spin_next_nbor_old
    E_total_new=-J1*sum_spin_nbor_new-J2*sum_spin_next_nbor_new       
 E_new=Ene+E_total_new-E_total_old 
 E_new = clamp(E_new, Energy[1], Energy[end])
#Eold_idx=min(E_number,max(Enr_old,Ene+4*rand(-4:4)))
#Enew_idx=min(E_number,max(Enr_old,E_new+4*rand(-4:4)))  
        Enr_old_idx = find_energy_index((Ene), Energy)
        Enr_new_idx = find_energy_index((E_new), Energy) 
       # Enr_old_idx = findfirst(==(Ene), Energy)
       # Enr_new_idx = findfirst(==(E_new), Energy) 
#print(S,'\n')
#println("Enr_old: $Ene, Enr_new: $E_new")
#println("Enr_old_idx: $Enr_old_idx, Enr_new_idx: $Enr_new_idx")
DS=S[Enr_old_idx]-S[Enr_new_idx]  

    if DS>0
       prob=1
    else
       prob=min(1,exp(DS))
    end
#(S)=exp^(ln S(DS))=S(DS) ; G(E)=exp^(S(E))
#println("Step $m: Enr_old=$Enr_old, Enr_new=$Enr_new, DS=$DS, prob=$prob")
        if rand() < prob                      # DS为正说明能量升高，为负说明能量降低   
        S[Enr_new_idx]      += f             # 更新态密度
        Hist[Enr_new_idx]   += 1             # 更新直方图
        spin[ns] = t_value;                  # 更新自旋
        Ene = E_new                          # 更新能量
        else
        S[Enr_old_idx]      += f             # 不更新 态密度
        Hist[Enr_old_idx]   += 1             # 不更新 直方图  
        end
        step += 1                            #判断 是否检查直方图平坦性 
        mmstep=nn  
        
      #  if step % mmstep == 0
       #     print(Hist,'\n')
       # end
        if (step) % mmstep == 0
            change= isflat(Hist,mmstep)
          #  print(change,'\n')
              if change==1
                f    *= 0.95 # 修改因子
                Hist .= 0    # 重置直方图 
                step  = 0
      #        print("mstep=",m,' ',"flat!,f=", exp(f),'\n')
              end
        end
        if f < 10^-8
        end
end  
    # mstep
    
   # S=S.-S[1].+log(3)
    
   return S, Hist, Energy,spin
end


function logsumexp(x)
    mx = maximum(x)
    mx + log(sum(exp.(x .- mx)))
end

function calculate_thermal_properties(ll::Int, Energy, S)
    nx = ll
    N=nx^2
    num_temps = 320
    parT_values = range(8.0,1.6,length=num_temps)
    E_ave = zeros(Float64, num_temps)
    E2_ave = zeros(Float64, num_temps)
    specific_heat = zeros(Float64, num_temps)
    
    for t in 1:num_temps
    #    T =0.2*t            # T/J1
        T =parT_values[t]             # T
        beta = 1.0 / T
        exponent = S .- beta .* Energy
        log_Z= logsumexp(exponent)                                    #放到对数空间处理配分函数
        E_avg = sum(Energy .* exp.(exponent .- log_Z))
        E2_avg = sum(Energy.^2 .* exp.(exponent .- log_Z))       
        E_ave[t] = E_avg
        E2_ave[t] = E2_avg
        specific_heat[t] = (E2_avg - E_avg^2) / (T^2 * N)
    end
    return E_ave/N, E2_ave, specific_heat
end

@time begin        
function main(ll::Int,J1,J2,mstep)
nx=ll;
nbor=latticesq(nx)
#spin=initspin_jitai(nx)
#spin=initspin(nx*nx)
#print("试验自旋：",spin,'\n')
spin=initspin(nx*nx)
S,Hist,Energy=Wanglandau(nx,J1,J2,nbor,spin,mstep)
E_avg,E2_avg, C_v_avg=calculate_thermal_properties(nx, Energy, S) 
return E_avg,E2_avg, C_v_avg,S,Hist,Energy
end
E_ave,E2_ave,C_v,S,Hist,Energy=main(8,-2,1,10000000)    
print("平均能量:",E_ave,'\n')  #根据配分函数计算能量 
#print(spin,'\n')
#print(E_MC,'\n')  #根据格点计算能量
#print(E2_MC,'\n')
print("态密度: ", S,'\n')
print("直方图: ", Hist,'\n')
print("能量区间: ", Energy,'\n')
end
