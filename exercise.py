
# 

# In[1]:


import numpy as np
import pandas as pd  # type: ignore
import matplotlib.pyplot as plt # type: ignore
import pymc as pm   # type: ignore
import arviz as az   # type: ignore


# ### Exercise 1 (max 3 points)
# 
# The file [pigs.csv](./pigs.csv) (Koen, Erin, Newton, Erica, & Ellington, E. Hance. (2022). Data for: Evaluating potential sources of invasive wild pigs in Ontario. https://doi.org/10.5061/dryad.5dv41ns6j) contains data about a population of pigs.
# 
# The columns have this meaning:
# 
#  - PigType - type of pig sighting (domesticated, wild boar, unknown)
#  - DETECTED - 0 - random location, 1 - pig sighting
#  - dist_boar - distance (meters) to nearest known premise with wild boar (meters)
#  - dist_pig - distance (meters) to nearest known premise with domestic pig (meters)
#  - borderMI - distance (meters) to border with Michigan (excluding the upper peninsula)
#  - borderNY - distance (meters) to border with New York (meters)
#  - borderQC - distance (meters) to border with Quebec
# 
# 
# Load the data in a pandas dataframe and make a `bool` column `pig_sighting` which is `True` iff the pig was detected from a pig sighting site. Use the first column as index.

# In[3]:

#########################################################################################################################################################
#replace 1 0 with True False
df = pd.read_csv("pigs.csv",)
df['DETECTED'] = df['DETECTED'].replace({1: True, 0: False})
df = df.drop(columns='Unnamed: 0')
df


#dtype
howell = pd.read_csv('https://raw.githubusercontent.com/rmcelreath/rethinking/master/data/Howell1.csv', sep=';', dtype={'male': bool})


df = pd.read_csv('Howell1.csv',sep=';')





with open ('sars.fasta','r') as f:
    ff= f.read()
    fff= ff.split('\n')
    fff.remove(fff[0])
    seq=''
    for i in fff:
        seq = seq + i



with open ('nc_051526_1.fasta')as f:
    ff = f.read()
    fasta= ff.split('\n')
    seq = ''
    fasta.remove(fasta[0])
    for i in fasta:
        seq = seq + ''.join(i)
    
seq



butterflies = pd.read_csv('butterfly_data.csv' )
butterflies





with open ('sars.fasta','r') as f:
    ff= f.read()
    fff= ff.split('\n')
    fff.remove(fff[0])
    seq=''
    for i in fff:
        seq = seq + i



df = pd.read_csv('data.csv',index_col=0)
df.head()





#########################################################################################################################



#dtype
howell = pd.read_csv('https://raw.githubusercontent.com/rmcelreath/rethinking/master/data/Howell1.csv', sep=';', dtype={'male': bool})

howell.head()

howell.describe(percentiles=[0.05, 0.95])

# ## Exercise 1
#
# Describe data for young people (less than 18 years old).

howell[howell.age < 18].describe()


#dropna

df_gauss=df.dropna(axis=0)


df[(df['PigType'] == 'wild boar') & (df['DETECTED'] == True)]['dist_boar'].mean()


df['bool']=df['DYRK1A_N']>0.45


df[(df['PigType'] == 'wild boar') & (df['DETECTED'] == True)]['dist_pig'].mean()


df[(df['PigType'].isin(['domesticated - pot bellied', 'domesticated'])) & (df['DETECTED'] == True)]['dist_boar'].mean()



#2 value in one column isin    


df[(df['Species'].isin(['Acacia tortilis', 'Acacia sp.','Acacia mellifera']))]['d13C'].mean()

df[(df['Species'].isin(['Acacia tortilis', 'Acacia sp.','Acacia mellifera']))]['d13C'].std()


df[(df['PigType'].isin(['domesticated - pot bellied', 'domesticated'])) & (df['DETECTED'] == True)]['dist_pig'].mean()

df['area'] = df.apply(lambda row: area(row['borderMI'], row['borderMI'], row['borderQC']), axis=1)
df

#ascending sampling
df_filtered = df[df['PigType'] == 'wild boar']
df_filtered = df_filtered.sort_values(by='dist_pig')
df_filtered = df_filtered.head(5)
print(df_filtered['borderMI'])





#datframe index
df = pd.DataFrame(possible_triplets(seq),index=possible_triplets(seq))
df.columns=['triplets']
df['occurrences']= df.apply(lambda x: triplet_occurance(seq, triplet=x['triplets']),axis=1)
df.drop('triplets', axis=1, inplace=True)
df.index.name = "triplets"
df.head()



df['occurance_is_even'] = df['occurrences'].apply(lambda x: x % 2 == 0)



#groupby

df2 = pd.DataFrame(df.groupby('age').length.mean())
df2['2'] =df.groupby('age').length.std()
df2.columns = ['mean_by_age','std_by_age']
df2


group = pd.DataFrame(butterflies.groupby("subarea"))

df.groupby(by='Behavior')['H3MeK4_N'].mean()

grouped = butterflies.groupby('year_month')

butterflies['avg_coll_dist'] = butterflies.groupby('year_month').apply(distance)


# In[38]:


milad = butterflies.groupby('year_month') 



df['gender']= df.apply(lambda x: gender_detect(x['dna']),axis=1)



df['a_twins']= df.apply(lambda x:repeat_detect(x['dna'],'A'),axis=1 )



col=butterflies["organic"].apply(bool)
butterflies.organic = col



# In[24]:


col2=butterflies["alternate_management"].apply(bool)
butterflies.alternate_management = col2



# In[25]:



butterflies['year_month'] = butterflies['year'].apply(str) + butterflies['month']


# In[32]:


butterflies['coord'] = list(zip(butterflies['x'], butterflies['y']))


# In[85]:






# In[62]:


butterflies.loc[0:,['year_month','coord']]


data["for_test"].all()==data['class'].all()



data['percentuals'] = data.occurances.apply(lambda x : ((x/1273)*100))


data = pd.DataFrame({"letter":countss, "occurances":values})
# In[91]:
list(data[data.percentuals > 5].letter)



data['percentuals'] = data.occurances.apply(lambda x : ((x/1273)*100))

data

data = pd.DataFrame({"letter":countss, "occurances":values})

data


male = df[df.male==1]
male.head()
tall_males = male[male.height>= 122.0]
female = df[df.male==0]
tall_females = female[female.height>= 122.0]
female.head()





df["w_dens"] = df.weight / df.height
df['expected_w_dens'] = df.apply(lambda x: w_dens_by_gender(x['age'], x['male']), axis=1)




even_mean= df[df['occurance_is_even']==True]['occurrences'].mean()
odd_mean=df[df['occurance_is_even']==False]['occurrences'].mean()




mean = df['occurrences'].mean()
std = df['occurrences'].std()
df['standardized'] = (df['occurrences'] - mean) / std
df.head()


# In[91]:


print(df['standardized'].mean())
print(df['standardized'].std())


df.filter(regex='ght') ----> will only keep thoes columns with 'ght' characters


############################################################################################################################################


fig, ax = plt.subplots(nrows=2,figsize =(6,8))
ax[0].hist(df[df['male']==True]['age'], bins='auto',density=True)
ax[0].hist(df[df['male']==False]['age'], bins='auto',density=True)

ax[1].hist(df[(df['male']==True) & (df['height'] > 1.2)]['age'], bins='auto',density=True)
ax[1].hist(df[(df['male']==False) & (df['height'] > 1.2)]['age'], bins='auto',density=True)



plt.hist(howell[(howell['male']==True) & (howell['age']>18)]['height']
         ,bins='auto',density=True,label='Adult Males')

plt.hist(howell[(howell['male']==False) & (howell['age']>18)]['height']
         ,bins='auto',density=True,label='Adult Femalse',alpha=.8)




#histogram
fig, ax = plt.subplots()
ax.hist(df.borderMI,bins= 'auto', density = True ,alpha = .75, label= 'Distance from Michigan', color='blue')
ax.hist(df.borderNY,bins= 'auto', density = True ,alpha = .75, label= 'Distance from NewYork', color='red')
ax.hist(df.borderQC,bins= 'auto', density = True ,alpha = .75, label= 'Distance from Quebec', color='Green')
plt.legend()


# In[6]:


#or histogram
fig = plt.figure(figsize = (24,7))
ax1 = fig.add_subplot(1, 3, 1)
ax1.hist(df.borderMI,bins= 'auto', density = True ,label= 'Distance from Michigan', color='blue')
ax2 = fig.add_subplot(1, 3, 2)
ax2.hist(df.borderNY,bins= 'auto', density = True ,label= 'Distance from NewYork', color='red')
ax3 = fig.add_subplot(1, 3, 3)
ax3.hist(df.borderQC,bins= 'auto', density = True ,label= 'Distance from Quebec', color='Green')
fig.legend()






fig, ax = plt.subplots()
plt.scatter(df[df['Plant_part']=='wood']['d13C'],df[df['Plant_part']=='wood']['d15N'],label='wood')
plt.scatter(df[df['Plant_part']=='leaves']['d13C'],df[df['Plant_part']=='leaves']['d15N'],label='leaves')
plt.scatter(df[df['Plant_part']=='bark']['d13C'],df[df['Plant_part']=='bark']['d15N'],label='bark')
plt.scatter(df[df['Plant_part']=='inflorescence']['d13C'],df[df['Plant_part']=='inflorescence']['d15N'],label='inflorescence')
plt.legend()
plt.show()


fig, ax = plt.subplots()
plt.scatter(df[df['DETECTED']==True]['dist_boar'],df[df['DETECTED']==True]['dist_pig'],color='darkred')
plt.scatter(df[df['DETECTED']==False]['dist_boar'],df[df['DETECTED']==False]['dist_pig'],color='gold')


#vertical mean line
fig, ax = plt.subplots(ncols=2,figsize =(10,5))
ax[0].hist(df[df['occurance_is_even']==True]['occurrences'], bins='auto',density=True,label = 'even triplets',color='darkred')
ax[0].axvline(even_mean, color='yellow', linestyle='--', label='Mean even')
ax[1].hist(df[df['occurance_is_even']==False]['occurrences'], bins='auto',density=True,label = 'odd triplets')
ax[1].axvline(odd_mean, color='orange', linestyle='--', label='Mean odd')
fig.legend()




fig, ax = plt.subplots(ncols=2,figsize=(10,5))
ax[0].hist(even_posterior['mu_even'], bins='auto', density=True, label='Posterior even means')
ax[0].hist(odd_posterior['mu_odd'], bins='auto', density=True, label='Posterior odd means',alpha=0.7)
ax[0].legend()
ax[1].hist(even_posterior['sigma_even'], bins='auto', density=True, label='Posterior sigma means')
ax[1].hist(odd_posterior['sigma_odd'], bins='auto', density=True, label='Posterior sigma means')
plt.legend()

# In[30]:


plt.hist(df.length, bins = 'auto', density=True,color = 'darkgreen' ,label='DNA Length ')
plt.legend()
plt.show()


plt.hist(df[df['gender']=='male']['length'],bins='auto',color='darkred',label='Length for male Sarchiapi E')
plt.legend()
plt.show()




plt.scatter(group[1][0].x, group[1][0].y )
plt.scatter(group[1][1].x, group[1][1].y )
plt.scatter(group[1][2].x, group[1][2].y )
plt.scatter(group[1][3].x, group[1][3].y )


plt.pie(countz(seq,1).values(),labels=countz(seq,1).keys())
plt.show()




plt.pie(countz(seq,1).values(),labels=countz(seq,1).keys())
plt.show()


plt.hist(df.age,bins='auto',density = False)
plt.show()


# In[7]:


#fig, ax = plt.subplots(ncols=2,figsize=(10,5))
fig = plt.figure(figsize = (10,5))
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)
ax1.hist(male.age,bins= 'auto', density = True ,label= 'males')
ax1.hist(female.age, bins = 'auto', density = True,label= 'females')
ax2.hist(tall_males.age,bins= 'auto', density = True, )
ax2.hist(tall_females.age, bins = 'auto', density = True)
fig.legend()




plt.scatter(df[df['male']==1]['w_dens'], df[df['male']==1]['age'], c='blue', label='Male')
plt.scatter(df[df['male']==0]['w_dens'], df[df['male']==0]['age'], c='red', label='Female')
plt.legend()





fig,ax=plt.subplots(nrows=8, figsize=(10,30))
for i in range(8):

  ax[i].scatter(data[data['class']==data["class"].unique()[i]]['standard_bacten'], data[data['class']==data["class"].unique()[i]]['standard_tau'], 
                label= data["class"].unique()[i],c='')
  ax[i].legend()





fig, ax =plt.subplots(nrows=78,ncols=2, figsize=(5, 3*78))
for i in range(1, 78):
  ax[i][0].hist(df[df['Genotype']=='Control'].iloc[:,i],density=True , bins='auto',label=df[df['Genotype']=='Control'].iloc[:,i].name)
  ax[i][0].legend()
  ax[i][0].hist(df[df['Genotype']=='Ts65Dn'].iloc[:,i],density=True , bins='auto',label=df[df['Genotype']=='Control'].iloc[:,i].name)
  ax[i][0].legend()
  ax[i][1].hist(df[df['Treatment']=='Memantine'].iloc[:,i],density=True , bins='auto',label=df[df['Genotype']=='Control'].iloc[:,i].name)
  ax[i][1].legend()
  ax[i][1].hist(df[df['Treatment']=='Saline'].iloc[:,i],density=True , bins='auto',label=df[df['Genotype']=='Control'].iloc[:,i].name)
  ax[i][1].legend()



def gaussian(x, mean, std):
    return (1 / (np.sqrt(2 * np.pi) * std)) * np.exp(-((x - mean) ** 2) / (2 * (std ** 2)))

# Load data into a Pandas DataFrame
df = pd.read_csv("data.csv")

# Loop through each row of the DataFrame
for index, row in df.iterrows():
    mean = row['mean']
    std = row['std']
    x = np.linspace(mean - 4 * std, mean + 4 * std, 100)
    y = gaussian(x, mean, std)
    
    # Create a new plot for each row
    fig, ax = plt.subplots()
    ax.plot(x, y)
    
    # Add labels and a title to the plot
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"Gaussian Function for row {index}")
    
    # Display the plot
    plt.show()


###########################################################################################################################################

model= pm.Model()
with model:
    alpha=pm.Normal('alpha', 0,5)
    beta=pm.Normal('beta',0,5)
    sigma=pm.Exponential('sigma', 1)
    mean=alpha + beta* df['C_cont']
    N_cont_obs=pm.Normal('N_count_obs',mean, sigma ,observed=df['N_cont'])
    trace=pm.sample()
    
pm.plot_posterior(trace)    

pm.summary(trace)    




borderQC_model = pm.Model()

with borderQC_model:
    
    mu = pm.Normal('mu_1',170000, 100000)
    sigma = pm.Exponential('sigma_1',1)
    borderQC = pm.Normal('borderQC',mu,sigma,observed=df['borderQC'])
    trace_1 = pm.sample()


# In[58]:


az.summary(trace_1)




#model creation and sampling posterior
stat_model_even = pm.Model()
with stat_model_even:
    mu_even = pm.Normal('mu_even',0,2)
    sigma_even = pm.Normal('sigma_even',0,2)
    even = pm.Normal('even',mu_even,sigma_even,observed=df[df['occurance_is_even']==True]['occurrences'])
    even_sample =pm.sample(chains=4)
even_posterior = az.extract(even_sample, combined=True).to_pandas()


# In[102]:


stat_model_odd = pm.Model()
with stat_model_odd:
    mu_odd = pm.Normal('mu_odd',0,2)
    sigma_odd= pm.Normal('sigma_odd',0,2)
    odd = pm.Normal('odd',mu_odd,sigma_odd,observed=df[df['occurance_is_even']==False]['occurrences'])
    odd_sample =pm.sample(chains=4)
odd_posterior = az.extract(odd_sample, combined=True).to_pandas()




model = pm.Model()

with model:
    mu = pm.Normal('mu',6.195)
    sigma = pm.Uniform('sigma', 0,10)
    hyp = pm.Normal('hypothes',mu,sigma)


# In[163]:


plt.hist(pm.draw(sigma,draws=10000),bins='auto')
plt.show()


# In[166]:


model2 = pm.Model()

with model2:
    mu2 = pm.Normal('mu_2',6.195)
    sigma2 = pm.Uniform('sigma_2', 0,10)
    hyp2 = pm.Normal('hypothes_2',mu2,sigma2,observed=df['a_twins'])

plt.hist(pm.draw(sigma2,draws=1000000),bins='auto')
plt.show()





model = pm.model()

with model:
    mu = pm.Normal('mu',0,5)
    
    sigma = pm.Exponential('std',1)
    
    



 model2 =pm.Model()
with model2:
  mu2=pm.Normal('mu2' , 2.15, 0.4)
  value2=pm.Normal('value2', mu2, observed=data['Bcatenin_N'])
     
plt.hist(pm.draw(value ,draws=10000), bins='auto')
plt.hist(pm.draw(value2 ,draws=10000), bins='auto')
plt.show()



model=pm.Model()
with model:
    mu=pm.Normal('mu', 0, 5)
    sigma=pm.Exponential('sigma', 1)
    observ=df['expected_w_dens'] - df['w_dens']
    error=pm.Normal('error', mu,sigma,observed=observ)
    trace=pm.sample()

pm.plot_posterior(trace)
pm.summary(trace)    
###################################################################################################################################################










howell[howell['age']>18]['height'].describe()


def Gaussian (data: np.ndarray, mu: float = 0 ,sigma: float = 1):
    return (1/(2*np.pi*sigma**2)**0.5)*np.exp(-((data-mu)**2)/(2*sigma**2))


x = np.linspace(-5,5,100)
fig , ax = plt.subplots()
_ = ax.plot(x, Gaussian(x))

# +

ax.plot(x,Gaussian(x,howell[howell['age']>18]['height'].mean()
                       ,howell[howell['age']>18]['height'].std()), label='Ideal')
_ = fig.legend()



fig , ax= plt.subplots(nrows=18, figsize=(6,20))
for i in range(18):
    ax[i].plot(x,Gaussian(x,df_gauss['d15N'].iloc[i],df_gauss['d15N_sd'].iloc[i]))




dicts = {}
keys = range(4)
values = ["Hi", "I", "am", "Milad"]
for i in keys:
        dicts[i] = values[i]
print(dicts)





def area(x:[float,int], y:[float,int], z:[float,int])-> float:
    """
    Computes the area of a circle with the biggest radius input
    
    >>> area(2000,3000,1000)
    28.274333882308138 
    """
    
    mac = max(x,y,z)
    mac = mac/1000
    area = np.pi*mac**2
    return area
area(2000,3000,1000)









def possible_triplets(sequence:str)-> list:
    i=0
    ls=[]
    while i < len(sequence):
        if len(sequence[i:i+3])<3:
            break
        ls.append(sequence[i:i+3])
        i +=1
    for i in ls:
        if i[::-1] in ls:
            ls.remove(i[::-1])
    set1 = set(ls)
    final = list(set1)
    return final
print(len(possible_triplets(seq)))
len(np.unique(possible_triplets(seq)))




def triplet_occurance (string:str, triplet:str)->int:
    """
    Returns the number of occurances of a triplet or it's reverse inside a string
    
    >>> string='CAATAATCC'
    >>> triplet='AAT'
    >>> triplet_occurance(string,triplet)
    3
    
    """
    cnt=0
    i=0
    while i < len(string):
        if len(string[i:i+3])<3:
                break
        if triplet == string[i:i+3] or triplet[::-1] == string[i:i+3] :
            cnt+=1
        i+=1
    return cnt
triplet_occurance(seq,'CAA')








def gender_detect(string:str)-> str:
    if string[0] == 'A' or string[0] == 'C':
        gender = 'male'
    else:
        gender = 'female'
    return gender
gender_detect('CTCCGTCTGTCCAGTACCTCTTC')



def repeat_detect(string:str, char:str)-> int:
    
    
    """Returns how many times the substring composed by the character repeated twice! 

    >>> string = 'ZXXZXXXZCCCX'
    >>> char = 'X'
    >>> repeat_detect(string,'char')
    3
    
    """
    cnt=0
    loop = 0
    while loop < len(string)-1:
        if string[loop] == char and string[loop+1]==char:
            cnt+=1
        loop+=1
    return cnt

repeat_detect(string,'X')












# In[29]:


def distance (x:int,x2:int,y:int,y2:int)-> float:
    """
    This function takes ...
    
    >>> Example:
    distance(2,3,5,6)
    1.4142135623730951
    """
    
    return math.sqrt((x-x2)**2 + (y-y2)**2)
distance(2,3,5,6)





def mk_class(a: str ,b: str, c:str) -> str:
  letter1=a[0].lower()
  b_string=''
  final=''
  for i in range(len(b)):
    if 65 < ord(b[i]) < 90 or 97 < ord(b[i]) < 123 :
      b_string=b_string + "".join(b[i])
  letter3=c[0].lower()
  final=letter1 + '-' + b_string + '-' + letter3 

  return final

























def euclidean_distance(x,x2,y,y2):
    return np.sqrt((x - x2)**2 + (y - y2)**2)


# In[93]:


for i in range(len(butterflies.coord)-1):
    x = butterflies.coord[i][0]
    y = butterflies.coord[i][1]
    x2 = butterflies.coord[i+1][0]
    y2 = butterflies.coord[i+1][1]
    ls.append(euclidean_distance(x,x2,y,y2))
    print (dis)



# In[ ]:


def avg_coll_dist ():
    lst = []
    
    
    return math.sqrt((x-x2)**2 + (y-y2)**2)









def countz(string: str, n :int) ->dict:
    
    lst=[]
    for i in range(0,len(string),n):
        lst.append(string[i:i+n])

    for i in lst:
        if len(i)%n !=0:
            lst.remove(i)


    dic={}
    for i in lst:
        dic[i]=lst.count(i)


    return dic




milad = countz(seq,2)
cnt = 0
if len(milad) == 276:
    cnt+=1
print (cnt)





milad_list = milad.tolist()


# In[153]:


milad2 = list(milad_list)





np.random.choice(milad2, size=1273)




np.random.seed(10)
ls = []
for i in range(10000):
    ls.append(''.join(np.random.choice(milad2, size=1273)))


# In[211]:


np.random.seed(10)
lss = []
for i in range(10000):
    lss.append(''.join(np.random.choice(milad2, size=1273)))


def countz(string: str, n :int) ->dict:
    
    lst=[]
    for i in range(0,len(string),n):
        lst.append(string[i:i+n])

    for i in lst:
        if len(i)%n !=0:
            lst.remove(i)


    dic={}
    for i in lst:
        dic[i]=lst.count(i)


    return dic

milad = countz(seq,2)
cnt = 0
if len(milad) == 276:
    cnt+=1
print (cnt)





milad2 = list(milad_list)

np.random.choice(milad2, size=1273)

np.random.seed(10)
ls = []
for i in range(10000):
    ls.append(''.join(np.random.choice(milad2, size=1273)))

np.random.seed(10)
lss = []
for i in range(10000):
    lss.append(''.join(np.random.choice(milad2, size=1273)))

### Exercise 7 (max 4 points)

Count how many strings among the ones produced during the previous exercise have at most the same number (276) of sub-sequences of 2 letters.

cnt=0
for i in lss:
    temp = countz(i,2)
    if len(temp) == 276:
        print ("TRUE")



def w_dens_by_gender (age, value:bool):
    if value == True:
        if age <= 30:
            res = age / 100
        elif age > 30 :
            res = age / 30
    elif value == False:
        if age <= 25:
            res = age / 110
        elif age > 25:
            res = age / 23
    return res
w_dens_by_gender(30,False)


# In[55]:

