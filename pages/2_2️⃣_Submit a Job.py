import torch
import streamlit as st
from transformers import (TrainingArguments,
                          Trainer,
                          GPT2Config,
                          GPT2Tokenizer,
                          AdamW, 
                          get_linear_schedule_with_warmup,
                          GPT2ForSequenceClassification)
import pickle
from Interpreter import Interpreter 
import numpy as np
import time
from stqdm import stqdm
import tkinter


st.write('### DDI-GPT: Predict DDIs using KG enhanced LLM')
# st.markdown("A demonstration using [Streamlit](https://streamlit.io/) with [HuggingFace's GPT-2](https://github.com/huggingface/transformers/).")


## input two drugs
drug1 = st.text_input('Drug1 eg: DB00343',value='DB00343')
drug2 = st.text_input('Drug2 eg: DB01268',value='DB01268')

uqcomb_dpwy_dict = pickle.load(open("./dict/uqcomb_dpwy_dict.p", "rb"))
uqcomb_dg_tg_dict = pickle.load(open("./dict/uqcomb_dg_tg_dict.p", "rb"))
uqcomb_dg_trans_dict = pickle.load(open("./dict/uqcomb_dg_trans_dict.p", "rb"))
uqcomb_dg_enzy_dict = pickle.load(open("./dict/uqcomb_dg_enzy_dict.p", "rb"))

# genes
sentence_dual = []
for drugname in [drug1, drug2]:
    sentence = []
    sentence.append(["Drug_is"])
    sentence.append([drugname])
    sentence.append(["pathway_is"])
    try:
        sentence.append(uqcomb_dpwy_dict[drugname][:2])
    except:
        sentence.append(";;")
    sentence.append(["target_is"])
    try:
        sentence.append(uqcomb_dg_tg_dict[drugname][:2])
    except:
        sentence.append(";;")
    sentence.append(["transporter_is"])
    try:
        sentence.append(uqcomb_dg_trans_dict[drugname][:2])
    except:
        sentence.append(";;")
    sentence.append(["enzyme_is"])
    try:
        sentence.append(uqcomb_dg_enzy_dict[drugname][:2])
    except:
        sentence.append(";;")
    
    flat_list = [item for sublist in sentence for item in sublist]
    sentence_dual.append(flat_list)

flat_duals = [item for sublist in sentence_dual for item in sublist]
text = ' '.join(flat_duals)

# text = st.text_input('Drug_is DB00788 pathway_is ; ; target_is HGNC:4932 transporter_is HGNC:42 HGNC:10959 enzyme_is HGNC:2625 HGNC:2637 Drug_is DB01095 pathway_is KEGG:map07027 KEGG:map07234 target_is HGNC:11059 HGNC:4624 transporter_is HGNC:40 enzyme_is HGNC:12530')
print(text)
st.markdown(text)
go = st.button('Predict')

model_name_or_path = 'cx229/pretrained'
n_labels = 2
model_config = GPT2Config.from_pretrained(pretrained_model_name_or_path=model_name_or_path, num_labels=n_labels)
# Get model's tokenizer.
tokenizer = GPT2Tokenizer.from_pretrained(pretrained_model_name_or_path="healx/gpt-2-pubmed-medium")
# default to left padding
tokenizer.padding_side = "left"
# Define PAD Token = EOS Token = 50256
tokenizer.pad_token = tokenizer.eos_token
# Get the actual model.
model = GPT2ForSequenceClassification.from_pretrained(pretrained_model_name_or_path=model_name_or_path, config=model_config)
# resize model embedding to match new tokenizer
model.resize_token_embeddings(len(tokenizer))
# fix model padding token id
model.config.pad_token_id = model.config.eos_token_id

if go and (drug1.startswith("DB")) and (drug2.startswith("DB")):
    try:
        max_sequence_len = None
        global input_ids
        input_ids = tokenizer.encode(text=text, return_tensors="pt", 
                                                padding=True, truncation=True,  
                                                max_length=max_sequence_len)

        with torch.no_grad(): 
            outputs = model(input_ids)
            logits, _  = outputs[:2]
            # Move logits and labels to CPU
            logits = logits.detach().cpu().numpy()
            # get predicitons to list
            predict_content = logits.argmax(axis=-1).flatten().tolist()
            if 0 in predict_content:
                result = "negative"
            else:
                result = "predicted interact"
            
            progress_text = "Calculation in progress. Please wait."
            for _ in stqdm(range(50), desc=progress_text):
                time.sleep(0.1)
                # my_bar.progress(percent_complete + 1, text=progress_text)

            st.markdown("GPT output: " + result)

    except Exception as e:
        st.exception("Exception: %s\n" % e)

elif (not drug1.startswith("DB")) | ( not drug2.startswith("DB")):
    st.markdown("input drug ids invalid")

fig = st.button('Measure token importance')
if fig:
        ## visualize
    def Phi(x):
        global model
        result = model(inputs_embeds=x)[0]
        return result # return the logit of last word
    input_embedding_weight_std = (
        model.get_input_embeddings().weight.view(1,-1)
        .std().item()
    )

    # maybe can global for two variables here
    max_sequence_len = None
    input_ids = tokenizer.encode(text=text, return_tensors="pt", 
                                    padding=True, truncation=True,  
                                    max_length=max_sequence_len)
    with torch.no_grad(): 
        x = model.get_input_embeddings()(input_ids).squeeze()

    revtext_list = []
    # text_dict = [revtext_list.append(tokenizer.decode(x)) for x in input_ids.tolist()[0]]
    for _ in input_ids.tolist()[0]:
        revtext_list.append(tokenizer.decode(_))

    interpreter = Interpreter(x=x, Phi=Phi, 
                            scale=10*input_embedding_weight_std,
                            words=revtext_list).to(model.device)
    with st.spinner('Wait for it...'):
        # This will take sometime.
        interpreter.optimize(iteration=5, lr=0.01, show_progress=False)
        sigma = interpreter.get_sigma()
    st.success('Done!')

    #https://stackoverflow.com/questions/62317723/tokens-to-words-mapping-in-the-tokenizer-decode-step-huggingface
    #start index because the number of special tokens is fixed for each model (but be aware of single sentence input and pairwise sentence input)
    tokens = tokenizer.tokenize(text)

    idx = 0
    enc =[tokenizer.encode(x, add_special_tokens=False, add_prefix_space=True) for x in text.split()]
    desired_output = []
    for token in enc:
        tokenoutput = []
        for ids in token:
            tokenoutput.append(idx)
            idx +=1
        desired_output.append(tokenoutput)

    highlevel_sigma = []
    for i in desired_output:
        sub_level = []
        for j in i:
            sub_level.append(sigma[j])
        highlevel_sigma.append(sum(sub_level))
            
    import matplotlib.pyplot as plt
    def visualize(sigma, words):
        """ Visualize the information loss of every word.
        """
        fi , ax = plt.subplots()
        ax.imshow([sigma], cmap='GnBu_r')
        ax.set_xticks(range(len(sigma)))
        ax.set_xticklabels(words)
        ax.set_yticks([0])
        ax.set_yticklabels([''])
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.show()
        st.pyplot(fi)

    visualize(sigma = highlevel_sigma, words = text.split())
    # st.set_option('deprecation.showPyplotGlobalUse', False)




# Footer
st.markdown(
    """
    <br>
    <h6>Chengqi Xu <a href="https://github.com/Mew233" target="_blank">GitHub Repo</a></h6>
    <h6>@ElementoLab, Weill Cornell Medicine</h6>
    """, unsafe_allow_html=True
    )