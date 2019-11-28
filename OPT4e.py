
# coding: utf-8

# In[2]:


#########################################OPT4e
###Optimal-features Predictor for T4SS Effector proteins
###by: Zhila Esna Ashari
    
import numpy as np
import pandas as pd
import os
import sys
import subprocess
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import scale
from sklearn.preprocessing import StandardScaler
import shutil
from shutil import copyfile
import tkinter
from tkinter import *
###################################################add a toggle button to show/hide the code##########################
from IPython.display import display
from IPython.display import HTML
import IPython.core.display as di 

di.display_html('<script>jQuery(function() {if (jQuery("body.notebook_app").length == 0) { jQuery(".input_area").toggle(); jQuery(".prompt").toggle();}});</script>', raw=True)

di.display_html('''<button onclick="jQuery('.input_area').toggle(); jQuery('.prompt').toggle();">Toggle code</button>''', raw=True)


####################################################------------------------------------------############################

def calculating_features(filename):
    
    
    # with open('seq3.fasta') as f:
    #     lines = f.readlines()
    
    fasta_file=filename + ".fa"
    lines = [line.rstrip('\n') for line in open(fasta_file)]

    # print(lines)

    
    my_sequence= ''
    sequences= []
    sequences_name= []

    # my_sequence_name=lines[0]

    for l in lines:
        if(l == ""): #ignore blank lines
                pass
        elif (l[0] == '>'):
            sequences_name.append(l[1:])
            sequences.append(my_sequence)
            my_sequence= ''
        else:
            my_sequence+= l

    sequences.append(my_sequence)
    del sequences[0]

    # sequences[1][0]
    # print(sequences)
    # print(sequences_name[0])
    # len(sequences_name)


    
    # output_features_all=np.zeros((len(sequences_name),1))
    output_features_all=pd.DataFrame()
    # output_features_all['names']=sequences_name
    # print(output_features_all)





    ####################################################------------------------------------------############################
    ###FUNCTIONS#######################
    

    def hydropathy(seq, num):
        hydropathy= 0

        if (num==1):
            seq2=seq
        elif(num==2):
            seq2=seq[-25:]
        elif (num==3):
            seq2=seq[0:25]

        for amino_acid in seq2:
            if (amino_acid=='G'):
                hydropathy-=0.4
            elif (amino_acid=='A'):
                hydropathy+=1.8
            elif (amino_acid=='P'):
                hydropathy-=1.6
            elif (amino_acid=='V'):
                hydropathy+=4.2
            elif (amino_acid=='L'):
                hydropathy+=3.8
            elif (amino_acid=='I'):
                hydropathy+=4.5
            elif (amino_acid=='M'):
                hydropathy+=1.9
            elif (amino_acid=='F'):
                hydropathy+=2.8
            elif (amino_acid=='Y'):
                hydropathy-=1.3
            elif (amino_acid=='W'):
                hydropathy-=0.9
            elif (amino_acid=='S'):
                hydropathy-=0.8
            elif (amino_acid=='T'):
                hydropathy-=0.7
            elif (amino_acid=='C'):
                hydropathy+=2.5
            elif (amino_acid=='N'):
                hydropathy-=3.5
            elif (amino_acid=='Q'):
                hydropathy-=3.5
            elif (amino_acid=='K'):
                hydropathy-=3.9
            elif (amino_acid=='H'):
                hydropathy-=3.2
            elif (amino_acid=='R'):
                hydropathy-=4.5
            elif (amino_acid=='D'):
                hydropathy-=3.5
            elif (amino_acid=='E'):
                hydropathy-=3.5
        return hydropathy

    ############################
    
    def Hydropathy_avg(seq, total_hydropathy):
        Hydropathy_avg= total_hydropathy/ len(seq)
        return Hydropathy_avg

    ############################

    
    def Polarity(seq):
        Polarity= 0.0

        for amino_acid in seq:
            if (amino_acid=='G'):
                Polarity+=9
            elif (amino_acid=='A'):
                Polarity+=8.1
            elif (amino_acid=='P'):
                Polarity+=8
            elif (amino_acid=='V'):
                Polarity+=5.9
            elif (amino_acid=='L'):
                Polarity+=4.9
            elif (amino_acid=='I'):
                Polarity+=5.2
            elif (amino_acid=='M'):
                Polarity+=5.7
            elif (amino_acid=='F'):
                Polarity+=5.2
            elif (amino_acid=='Y'):
                Polarity+=6.2
            elif (amino_acid=='W'):
                Polarity+=5.4
            elif (amino_acid=='S'):
                Polarity+=9.2
            elif (amino_acid=='T'):
                Polarity+=8.6
            elif (amino_acid=='C'):
                Polarity+=5.5
            elif (amino_acid=='N'):
                Polarity+=11.6
            elif (amino_acid=='Q'):
                Polarity+=10.5
            elif (amino_acid=='K'):
                Polarity+=11.3
            elif (amino_acid=='H'):
                Polarity+=10.4
            elif (amino_acid=='R'):
                Polarity+=10.5
            elif (amino_acid=='D'):
                Polarity+=13
            elif (amino_acid=='E'):
                Polarity+=12.3
        return Polarity



    ############################


    
    def Mr(seq):
        Mr= 0


        for amino_acid in seq:
            if (amino_acid=='G'):
                Mr+=+75.07
            elif (amino_acid=='A'):
                Mr+=+89.09
            elif (amino_acid=='P'):
                Mr+=+115.13
            elif (amino_acid=='V'):
                Mr+=+117.15
            elif (amino_acid=='L'):
                Mr+=+131.17
            elif (amino_acid=='I'):
                Mr+=+131.17
            elif (amino_acid=='M'):
                Mr+=+149.21
            elif (amino_acid=='F'):
                Mr+=+165.19
            elif (amino_acid=='Y'):
                Mr+=+181.19
            elif (amino_acid=='W'):
                Mr+=+204.24
            elif (amino_acid=='S'):
                Mr+=+105.09
            elif (amino_acid=='T'):
                Mr+=+119.12
            elif (amino_acid=='C'):
                Mr+=+121.15
            elif (amino_acid=='N'):
                Mr+=+132.12
            elif (amino_acid=='Q'):
                Mr+=+146.15
            elif (amino_acid=='K'):
                Mr+=+146.19
            elif (amino_acid=='H'):
                Mr+=+155.16
            elif (amino_acid=='R'):
                Mr+=+174.20
            elif (amino_acid=='D'):
                Mr+=+133.10
            elif (amino_acid=='E'):
                Mr+=+147.13
        return Mr


    ##########################################################
    
    def aa_composition(seq,aa_list):
        aa_comp= [None] * 20
        for i in range(20):
            aa_comp[i]=seq.count(aa_list[i])/len(seq)

        return aa_comp
    #########################################################
    
    def dp_composition(seq,dp_list):
        dp_comp= [None] * 400
        for i in range(400):
            dp_comp[i]=seq.count(dp_list[i])/(len(seq)-1)
        return dp_comp


    ############################################################
    

    def pssm_composition(seq,aa_list,mypssm):
        pssm_comp= [0] * 400
    #     print(len(pssm_comp))
        for i in range(20):
            mysum=[0]*20
            char_index=[]
            for j in range(len(seq)):
                if aa_list[i]==seq[j]:
                    char_index.append(j)

            for k in range(len(char_index)):
                mysum=np.add(mysum,list(map(int, mypssm[char_index[k]])))

            pssm_comp[20*i:20*(i+1)-1]=mysum
    #         print(len(pssm_comp))
    #         print(mysum)
    #         print(pssm_comp)

        pssm_comp=pssm_comp[0:400]
        pssm_comp=np.divide(pssm_comp,float(len(seq)))

        return pssm_comp


    ############################################################
    

    def pssm_autocov(seq,aa_list,mypssm):
        G=8
        pssm_auto= [0] * 20 *G  ##########################G=8

        for i in range(20):
            t=[0]*G
            smean=sum(list(map(int, mypssm[:,i])))/len(seq)

            for j in range(G):
                for k in range(len(seq)-j-1):
                    J=(int(mypssm[k,i])-smean)*(int(mypssm[k+j+1,i])-smean)
                    t[j]=t[j]+J
                t[j]=t[j]/(len(seq)-j)
    #         print(t)
            pssm_auto[8*i:8*(i+1)-1]=t



        pssm_auto=pssm_auto[0:160]
        return pssm_auto

    ############################################################
    ####################################################------------------------------------------############################
    #############################Build up output columns############################
    for i in range(1,401):
        output_features_all["V%s"%i]=""

    for i in range(1,161):
        output_features_all["W%s"%i]=""

    # print(output_features_all)



    
    aa_list='ARNDCQEGHILKMFPSTWYV'
    for i in aa_list:
        output_features_all["%s"%i]=""

    dp_list= [None] * 400
    n=0
    for i in range(20):
        for j in range(20):
            a=aa_list[i]+aa_list[j]
            dp_list[n]=a
            n+=1
            output_features_all["%s"%a]=""

    # print(dp_list)
    # print(len(dp_list))

    output_features_all["total_hydropathy"]=""
    output_features_all["hydropathy_C_terminal"]=""
    output_features_all["hydropathy_N_terminal"]=""
    output_features_all["hydropathy_average"]=""
    output_features_all["polarity"]=""
    output_features_all["molecular_mass"]=""
    output_features_all["homology"]=""
    ##############run functions###################################################################

    seqs_length=[]
    total_hydropathy=[]
    hydropathy_C_terminal=[]
    hydropathy_N_terminal=[]
    hydropathy_average=[]
    polarity=[]
    molecular_mass=[]
    homology=[]


    ##################Run and loop over sequences#########################################################################
    counter=0
    counting_seq=0

    for seq in sequences:

    ################################################
        seqs_length.append(len(seq))

        total_hydropathy.append(hydropathy(seq, 1))
        hydropathy_C_terminal.append(hydropathy(seq, 2))
        hydropathy_N_terminal.append(hydropathy(seq, 3))
        hydropathy_average.append(Hydropathy_avg(seq, total_hydropathy[counter]))
        polarity.append(Polarity(seq))
        molecular_mass.append(Mr(seq))
        counter+=1

    ##########################################

        aa_comp=aa_composition(seq,aa_list)
        dp_comp=dp_composition(seq,dp_list)
        output_features_all.at[counting_seq,"A":"V"]=aa_comp
        output_features_all.at[counting_seq,"AA":"VV"]=dp_comp

    #     print(aa_comp)
    # #     print(dp_comp)
    #     print(seq)
    #     print(len(seq))

    ###***PSSM***###
    ################################################
    #####################################################################################################
    
        x_len=len(seq)
        ###run psi-blast to get pssm file
        myfile = open("seq_pssm.fa", "w")
        myfile.write(seq)
        myfile.close()
        subprocess.call('.\\blast\\bin\\blastpgp -i .\\seq_pssm.fa -d uniprot -j 3 -h 0.001 -a 8 -Q .\\outname_test.pssm', shell=True)

        ####read and parse and extract pssm 
        xx=open("outname_test.pssm", "r")
        yy=[line.rstrip('\n') for line in xx]  

        thefile = open('mytemp.txt', 'w')
        for item in yy[3:-6]:
            thefile.write( item )

        # thefile = open('mytemp.txt', 'w')#//////////////////////\\\\\\\\\\\\\\\\\
        # for item in yy[3:-6]:
        #     thefile.write( item )

        thefile.close()#//////////////////\\\\\\\\\\\\\\\\\\\


        input1 = np.loadtxt('mytemp.txt', dtype=(str, float ))
    #     print(input1)
    #     print(len(input1))
        myarray=np.array(input1)
        mypssm1=np.resize(myarray,(x_len, 44))
        mypssm=mypssm1[:,2:22 ]

    #     print(mypssm)
#         print(len(mypssm))
#         print(input1[-43])
        ################################create pssm-comp
        pssm_comp=pssm_composition(seq,aa_list,mypssm)
    #     print(pssm_comp)
    #     print(len(pssm_comp))
        ##############################create pssm-auto
        pssm_auto=pssm_autocov(seq,aa_list,mypssm)
    #     print(pssm_auto)
    #     print(len(pssm_auto))

    ######################
        output_features_all.at[counting_seq,"V1":"V400"]=pssm_comp
        output_features_all.at[counting_seq,"W1":"W160"]=pssm_auto

        counting_seq+=1


    ################################################################
        myfile = open("seq_hom.fa", "w")
        myfile.write(seq)
        myfile.close()

        homol=0
        subprocess.call('.\\ncbi-blast-2.7.1+\\bin\\blastp -query .\\seq_hom.fa -db effector_db.txt -out .\\outname.temp -evalue 0.01', shell=True)

        if '>' in open('outname.temp').read():
            homol=1
#         print (homol)
        homology.append(homol)
    
    ################################################################
        subprocess.call("pepcoil -sequence .\\input_sequences.fa -window 28 -coil -frame -noother -rformat motif -auto -outfile input_emboss.pepcoil")
        ##Coiled Coils
        #####Gather coiled coil regions from output file
        coiled_coils=[]
        with open('input_emboss.pepcoil') as f:###outfile is the output of pepcoil software
            lines1 = f.readlines()
        ####################################       
        for mylines in lines1:
            if "HitCount: " in mylines:
                coiled_coils.append(mylines[12])


        #         print(mylines[12])
        #         counter+=1
        # print(coiled_coils)

        # output_features_all["Coiled_coils"]=coiled_coils
        # print(output_features_all)

    #     pd.DataFrame(coiled_coils).to_csv("Coiled_coils_Anap_HGE_merged.csv", header=False, index=False)

        #######################################################
        signalp=[]
        with open('sig_out') as f2:###outfile is the output of pepcoil software
            lines2 = f2.readlines()
        ############################################################the specific number we are looking for is in the index 12 of the line       
        for mylines in lines2:
        #     print(mylines)
            if " N " in mylines:
                signalp.append(0)
            elif " Y " in mylines:
                signalp.append(1)

        # # print(signalp)
        # pd.DataFrame(signalp).to_csv("signalp_train3.csv", header=False, index=False)
        # # output_features_all["SignalP"]=signalp
        # # print(output_features_all)


    ####################################################------------------------------------------############################
    ################################ Edit output (add or remove features)
    # print(seqs_length)

    output_features_all["length"]=seqs_length
    # print(output_features_all)


    # print (total_hydropathy)
    # print (hydropathy_C_terminal)
    # print (hydropathy_N_terminal)
    # print (hydropathy_average)
    # print (polarity)
    # print(molecular_mass)


    output_features_all["total_hydropathy"]=total_hydropathy
    output_features_all["hydropathy_C_terminal"]=hydropathy_C_terminal
    output_features_all["hydropathy_N_terminal"]=hydropathy_N_terminal
    output_features_all["hydropathy_average"]=hydropathy_average
    output_features_all["polarity"]=polarity
    output_features_all["molecular_mass"]=molecular_mass
    output_features_all["homology"]=homology
    output_features_all["coiled_coil"]=coiled_coils
    # print(output_features_all)



    ########################
    #######

    ###create a list of them to remove (this is for dp)
    removing= []
    selected=["KE","EK","AG","KK","GG","VG","EE","KD","DK","KL","GL","VI","VV","GA","NN","GV","FF","LP","AV","VR","KI","AR","NS","IA","LA","WQ","PL","LE","RA","AW","EL","FS","LG","DS","SS","QV","SD","NL","TN","GS","IK","VW","AA","MA","MG","YL"]

    n=0
    for i in range(20):
        for j in range(20):
            a=aa_list[i]+aa_list[j]
            if a not in selected:
                removing.append(a)
            n+=1

    ######

    output_features_all=output_features_all.drop(['C', 'Q','H','I','L','P','T','W','Y'], axis=1) ###remove for aa
    output_features_all=output_features_all.drop(removing, axis=1)###remove for dp

    # print(output_features_all.shape)   
    # print(output_features_all)

    ########################

    ########################
    ###
    removing_pssm_comp= []
    removing_pssm_auto= []
    selected_pssm_comp=["V5","V303","V246","V113","V143","V49","V56","V15","V205","V397","V20","V232","V158","V323","V17","V23","V260","V156","V369","V308","V227","V278","V132","V399","V43","V144","V126","V136","V312","V41","V179","V2","V340","V222","V8","V29","V385","V12","V226","V236","V306","V251","V76","V145","V230","V304","V221","V10","V7","V33","V6","V155","V395","V124","V381","V237","V386","V228","V13","V279","V16","V392","V127","V149","V229","V224","V378","V393","V307","V387","V14","V154","V253","V18","V11","V57","V159","V142","V42","V280","V384","V234","V48","V19","V147","V4","V9","V240","V121","V67","V231","V26","V160","V286","V239","V27","V388","V250","V274","V60","V146","V225","V271","V379","V141","V153","V223","V383","V64","V238","V374","V107","V44","V63","V389","V66","V51","V58","V125","V152","V157","V398","V31","V37","V62","V150","V61","V50","V382","V297","V35","V131","V151","V130","V194","V3","V55","V248","V28","V396","V254","V321","V128","V54","V38","V25","V275","V137","V52","V112","V261","V242","V72","V241","V285","V68","V269","V252","V284","V21","V348","V293","V258","V39","V287","V266","V244","V292","V346","V361","V45","V140","V262","V34","V134","V288","V277","V371","V366","V394","V289","V71","V53","V133","V245","V59","V300","V257","V249","V214","V377","V36","V75","V391","V350","V148","V129","V122","V317","V301","V316","V69","V165","V400","V178","V256","V175","V65","V174","V77","V70","V243","V168","V259","V167","V247","V74","V79","V73","V233","V80","V78","V200","V185","V196","V195","V183","V184","V192","V188","V182","V198","V181"]
    selected_pssm_auto=["W12","W92","W113","W11","W91","W51","W44","W43","W52","W49","W25","W89","W41","W27","W68","W115","W10","W114","W9","W140","W66","W90","W57","W35","W67","W137","W26","W142","W42","W139","W18","W4","W65","W138","W28","W34","W116","W19","W50","W17","W1","W122","W20","W148","W33","W13","W124","W130","W121","W93","W36","W141","W58","W108","W88","W53","W45","W145","W132","W32","W37","W157","W129","W24","W104"]
    n=0
    for i in range(1,401):
        a="V%s"%i
        if a not in selected_pssm_comp:
            removing_pssm_comp.append(a)
    # print(removing_pssm_comp)

    for j in range(1,161):
        a="W%s"%j
        if a not in selected_pssm_auto:
            removing_pssm_auto.append(a)
    # print(removing_pssm_auto)
    ######

    output_features_all=output_features_all.drop(removing_pssm_comp, axis=1)###remove for pssm_comp
    output_features_all=output_features_all.drop(removing_pssm_auto, axis=1)###remove for pssm_auto

#     print(output_features_all.shape)   
#     print(output_features_all)
    output_features_all.to_csv("%s.csv" %filename, header=True, index=False)
    
    
#########################################     Prediction     ######################################  

def predict_candidate_effectors():
    
    ###READ trainig and test set and extract lables
    training_set1=pd.read_csv('training_set.csv')
    training_set=training_set1.iloc[:, 0:364]
    label=training_set1.iloc[:,364]
    training_rows=training_set.shape[0]
    
    testing_set=pd.read_csv('input_sequences.csv')
    testing_rows1=testing_set.shape[0]
    
    if testing_rows1<21:
        for i in range(training_rows):
            testing_set=testing_set.append(training_set.iloc[i,:])
        testing_set=testing_set.reset_index(drop=True)


    ##NORMALIZE
    scaler = StandardScaler()
    training_set = pd.DataFrame(scaler.fit_transform(training_set), index=training_set.index, columns=training_set.columns)
    testing_set = pd.DataFrame(scaler.transform(testing_set), index=testing_set.index, columns=testing_set.columns)
    
    
    if testing_rows1<21:
        testing_rows=testing_set.shape[0]
        for i in range(testing_rows1,testing_rows):
            testing_set=testing_set.drop(i)

    ###Predict effectors using machine learning

    svm_model_linear = SVC(kernel = 'rbf', C = 1).fit(training_set, label)
    predictions = svm_model_linear.predict(testing_set)

    svm_model_linear = SVC(kernel = 'linear', C = 1).fit(training_set, label)
    predictions2 = svm_model_linear.predict(testing_set)

    logistic_model=LogisticRegression().fit(training_set, label)
    predictions3 = logistic_model.predict(testing_set)

#     pd.DataFrame(predictions).to_csv("predict_svmrbf.csv", header=True, index=False)
#     pd.DataFrame(predictions2).to_csv("predict_svmlinear.csv", header=True, index=False)
#     pd.DataFrame(predictions3).to_csv("predict_LR.csv", header=True, index=False)

#########################################################################################
#########################################################################################
#########################################################################################
#################################################################
    lines = [line.rstrip('\n') for line in open('training_set_non-effectors.fasta')]
    # lines = [line.rstrip('\n') for line in open('Anaplasma_Phagocytophilum_effectors.fasta')]
    # print(lines)

    ###read all sequences ##############################################
    my_sequence= ''
    sequences_noneff= []
    sequences_name_noneff= []


    for l in lines:
        if(l == ""): #ignore blank lines
                pass
        elif (l[0] == '>'):
            sequences_name_noneff.append(l[1:])
            sequences_noneff.append(my_sequence)
            my_sequence= ''
        else:
            my_sequence+= l

    sequences_noneff.append(my_sequence)
    del sequences_noneff[0]


# #################################################################
    lines = [line.rstrip('\n') for line in open('input_sequences.fa')]
    # lines = [line.rstrip('\n') for line in open('Anaplasma_Phagocytophilum_effectors.fasta')]
    # print(lines)

    ###read all sequences in sequences and their names in sequences_name##############################################
    my_sequence= ''
    sequences_input= []
    sequences_name_input= []


    for l in lines:
        if(l == ""): #ignore blank lines
                pass
        elif (l[0] == '>'):
            sequences_name_input.append(l[1:])
            sequences_input.append(my_sequence)
            my_sequence= ''
        else:
            my_sequence+= l

    sequences_input.append(my_sequence)
    del sequences_input[0]





#########################################################################################
#########################################################################################
#########################################################################################
##################################################### (predictions analysis)#################################
    predict_file = open("predicted_effectors_total.fasta", "w+")

    count_effector=0   
    for i in range(len(predictions)):
        if (predictions[i]==1 and sequences_input[i] not in sequences_noneff):
            predict_file.write(str(sequences_name_input[i])+ '\n'+ '\n')
            count_effector+=1
    predict_file.close()
    del predict_file
#     print(count_effector)

    one_line="number of predicted effectors(total)="+str(count_effector) + "\n"  +"\n" 
    with open("predicted_effectors_total.fasta", 'r+') as fp:
        lines = fp.readlines()     # lines is list of line, each element '...\n'
        lines.insert(0, one_line)  # you can use any index if you know the line index
        fp.seek(0)                 # file pointer locates at the beginning to write the whole file again
        fp.writelines(lines)  


    #########################  see which one has most probability based on voting
    predict_file1 = open("predicted_effectors_mostlikely.fasta", "w+")
    predict_file2 = open("predicted_effectors_probable.fasta", "w+")
    predict_file3 = open("predicted_effectors_lesslikely.fasta", "w+")

    count_effector_most=0
    for i in range(len(predictions)):
        if (predictions[i]==1 and sequences_input[i] not in sequences_noneff):
            if (predictions2[i] ==1 and predictions3[i]==1):
                predict_file1.write(str(sequences_name_input[i])+ '\n'+ '\n')
                count_effector_most+=1


    predict_file1.close()
    del predict_file1
#     print("most:")
#     print(count_effector_most)

    one_line="number of predicted effectors(most likely)="+str(count_effector_most) + "\n"  +"\n" 
    with open("predicted_effectors_mostlikely.fasta", 'r+') as fp:
        lines = fp.readlines()     # lines is list of line, each element '...\n'
        lines.insert(0, one_line)  # you can use any index if you know the line index
        fp.seek(0)                 # file pointer locates at the beginning to write the whole file again
        fp.writelines(lines) 
    ###

    count_effector_mid=0
    for i in range(len(predictions)):
        if (predictions[i]==1 and sequences_input[i] not in sequences_noneff):
            if ((predictions2[i] ==1 and predictions3[i]==-1) or (predictions2[i] ==-1 and predictions3[i]==1)):
                predict_file2.write(str(sequences_name_input[i])+ '\n'+ '\n')
                count_effector_mid+=1

    predict_file2.close()
    del predict_file2
#     print("mid:")
#     print(count_effector_mid)
    one_line="number of predicted effectors(probable)="+str(count_effector_mid) + "\n"  +"\n" 
    with open("predicted_effectors_probable.fasta", 'r+') as fp:
        lines = fp.readlines()     # lines is list of line, each element '...\n'
        lines.insert(0, one_line)  # you can use any index if you know the line index
        fp.seek(0)                 # file pointer locates at the beginning to write the whole file again
        fp.writelines(lines) 
    ###

    count_effector_less=0
    for i in range(len(predictions)):
        if (predictions[i]==1 and sequences_input[i] not in sequences_noneff):
            if (predictions2[i] ==-1 and predictions3[i]==-1):
                predict_file3.write(str(sequences_name_input[i])+ '\n'+ '\n')
                count_effector_less+=1

    predict_file3.close()
    del predict_file3
#     print("less:")
#     print(count_effector_less)  
    one_line="number of predicted effectors(less likely)="+str(count_effector_less) + "\n"  +"\n" 
    with open("predicted_effectors_lesslikely.fasta", 'r+') as fp:
        lines = fp.readlines()     # lines is list of line, each element '...\n'
        lines.insert(0, one_line)  # you can use any index if you know the line index
        fp.seek(0)                 # file pointer locates at the beginning to write the whole file again
        fp.writelines(lines) 


#########################################     GUI     ######################################

top = tkinter.Tk()
top.title("OPT4e")
top.geometry("500x470")
top.configure(background="blue")


canvas = Canvas(top, width=235, height = 139)      
canvas.pack(padx=5, pady=5)      
# # photo1=PhotoImage(file=)
img = PhotoImage(file="logo.gif")      
canvas.create_image(120,71, anchor=CENTER, image=img) 

##############################

# # Label(top, image=photo1, bg="black").grid(row=0, column=0, sticky=w)
# thelable=Label(top, text="OPT4e", bg="royalblue", fg="white", font="none 20 bold")
# # .grid(row=0,column=0)
# thelable.pack()
# Code to add widgets will go here...

# window = Toplevel(top)###open a new window when a button is pressed

# # root = Tk()
T = Text(top,bg="dark blue", fg="white", font="none 12 bold",bd=7,height=2, width=33, padx=10, pady=10)
# insert_centered(T, 'Thinking...')


T.pack(padx=10, pady=10)
T.insert(END, "Optimal-features Predictor for T4SS                         Effector proteins   ","center")
# # # mainloop()


def create_window1():
    window = Toplevel(top)
    window.geometry("320x60")
    window.configure(background="grey")
    L1 = Label(window, bg="grey",fg="white",text="Results stored in predicted_effectors... files \nThank you")
    L1.pack( )
    
    B1=Button(window, text="Okay", command=window.destroy) 
    B1.pack()
    
    
def create_window2():
    window = Toplevel(top)
    window.geometry("320x60")
    window.configure(background="grey")
    L1 = Label(window, bg="grey",fg="white",text="New effector sequences added to software \nThank you")
    L1.pack( )
    
    B1=Button(window, text="Okay", command=window.destroy) 
    B1.pack()
    
def create_window3():
    window = Toplevel(top)
    window.geometry("320x60")
    window.configure(background="grey")
    L1 = Label(window, bg="grey",fg="white",text="New non-effector sequences added to software \nThank you")
    L1.pack( )
    
    B1=Button(window, text="Okay", command=window.destroy) 
    B1.pack()
    
def create_window4():
    window = Toplevel(top)
    window.geometry("320x60")
    window.configure(background="grey")
    L1 = Label(window, bg="grey",fg="white",text="Training set restored to original version. \nThank you")
    L1.pack( )
    
    B1=Button(window, text="Okay", command=window.destroy) 
    B1.pack()
    
    
    
    

def first_button():
    calculating_features('input_sequences')
    predict_candidate_effectors()
    create_window1()
  
    
def second_button():
    
    file_new = open("input_sequences.fa", "r")
    data_new = file_new.read()
    file_new.close()

    file_all= open("training_set_effectors.fasta","a")
    file_all.write("\n")
    file_all.write(data_new)
    file_all.write("\n")
    file_all.close()
    
    file_all_db= open("effector_db.txt","a")
    file_all_db.write("\n")
    file_all_db.write(data_new)
    file_all_db.write("\n")
    file_all_db.close()
    
    calculating_features('input_sequences')
    
    
    df1 = pd.read_csv('input_sequences.csv')
    l = [1] * df1.shape[0]
    df1["label"]=l
    df1.to_csv("seq_NoHeader.csv", header=False, index=False)
    
    file_new_features = open("seq_NoHeader.csv", "r")
    data_new_features = file_new_features.read()
    file_new_features.close()
    
    file_all_features= open("training_set.csv","a")
    file_all_features.write(data_new_features)
    file_all_features.close()
    ###ONE TIME####
    ###build up database for BLAST and homology
    subprocess.call('.\\ncbi-blast-2.7.1+\\bin\\makeblastdb -in effector_db.txt -dbtype prot', shell=True)
    create_window2() 
    
    
    
def third_button():
    file_new = open("input_sequences.fa", "r")
    data_new = file_new.read()
    file_new.close()
    
    file_all= open("training_set_non-effectors.fasta","a")
    file_all.write("\n")
    file_all.write(data_new)
    file_all.write("\n")
    file_all.close()
    
    calculating_features('seq')
    
    
    df1 = pd.read_csv('seq.csv')
    l = [-1] * df1.shape[0]
    df1["label"]=l
    df1.to_csv("seq_NoHeader.csv", header=False, index=False)
    
    file_new_features = open("seq_NoHeader.csv", "r")
    data_new_features = file_new_features.read()
    file_new_features.close()
    
    file_all_features= open("training_set.csv","a")
    file_all_features.write(data_new_features)
    file_all_features.close()
        
    create_window2() 
    
    
def fourth_button():
#     os.remove("seq3.csv")
#     os.remove("seq3.fasta")
    copyfile("training_set_effectors_original.fasta", "training_set_effectors.fasta")
    copyfile("training_set_non-effectors_original.fasta", "training_set_non-effectors.fasta")
    copyfile("training_set_original.csv", "training_set.csv")
    copyfile("effector_db_original.txt", "effector_db.txt")
    subprocess.call('.\\ncbi-blast-2.7.1+\\bin\\makeblastdb -in effector_db.txt -dbtype prot', shell=True)
    create_window4()  
    
# b = Button(text="Create new window", command=create_window)
# b.pack()

w = tkinter.Button ( text ="Predict effector candidates \n(provide sequences in input_sequences.fa)", activebackground="grey",command=first_button )
# .grid(row=1, column=1, padx=10, pady=10)
# w .grid(row=0, column=1, padx=10, pady=10)
w.pack(padx=10, pady=10)

w2 = tkinter.Button ( text ="Add new validated effectors to software \n(provide sequences in input_sequences.fa)" , activebackground="grey",command=second_button)
# .grid(row=2, column=1, padx=10, pady=10)
w2.pack(padx=10, pady=10)

w3 = tkinter.Button ( text ="Add new known non-effectors to software \n(provide sequences in input_sequences.fa)" , activebackground="grey",command=third_button)
# .grid(row=2, column=1, padx=10, pady=10)
w3.pack(padx=10, pady=10)

w4 = tkinter.Button ( text ="Reset software to original version", activebackground="grey",command=fourth_button)
w4.pack(padx=10, pady=10)





top.mainloop()


