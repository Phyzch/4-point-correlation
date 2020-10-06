//
// Created by phyzch on 6/20/20.
// This file include function to compute system reduced density matrix fast.
//

#include<iostream>
#include"system.h"
#include"util.h"
using namespace std;

int compare_vmode (const vector<int> & vmode1, const  vector<int> & vmode2, const vector<int> & vmode1_0, const vector<int> & vmode2_0){
    // first compare vmode1, in this way, we don't have to merge sort between different process.
    int size1= vmode1.size();
    if( size1 != vmode1_0.size()){
        cout<< "error! two vector's mode should equal."<<endl;
        exit(-8);
    }
    int size2= vmode2.size();
    if(size2 != vmode2_0.size()){
        cout<< "error! two vector's mode should equal."<<endl;
        exit(-8);
    }

    int i;
    for(i=0;i<size1;i++){
        if(vmode1[i]>vmode1_0[i]) return 1;
        else if (vmode1[i] < vmode1_0[i]) return -1;
    }
    for(i=0;i<size2;i++){
        if(vmode2[i]>vmode2_0[i]) return 1;
        else if (vmode2[i]<vmode2_0[i]) return -1;
     }
    return 0;
}
// -------------- use binary esarch insert method to construct sys_state. This method is relatively slow. ----------------------
int find_location_binarysearch_sys_quotient_state(const vector<sys_quotient_state> & list, vector<int> & vmode1, const vector  <int> & vmode2, bool & exist){
    // used to sort quotient state list. return location to insert or manipulate, also return exist to indicate this quotient_state exist in list or not
    // use information of sys_state and vmode to locate the state.
    // use binary search method.
    int right_flag= list.size()-1; // this step is crucial for binary search.
    int left_flag= 0;
    if(list.size()==0){
        exist= false;
        return 0;
    }
    int position= (left_flag + right_flag)/2;
    exist =false;
    int mark;
    while(left_flag < right_flag){
        mark = compare_vmode(vmode1,vmode2,list[position].vmode1,list[position].vmode2);
        if(mark== 1){
            // state > list[position]
            left_flag=position+1;
            position= (left_flag + right_flag)/2;
        }
        else if (mark == -1){
            right_flag= position -1;
            position = (left_flag + right_flag )/2;
        }
        else if(mark==0){
            exist= true;
            return position;
        }
    }
    // Now situation is left_flag == right_flag == position or position == right_flag == left_flag -1. either case we compare list[position] with our state
    mark = compare_vmode(vmode1,vmode2,list[position].vmode1,list[position].vmode2);
    if( mark == 0 ){
        exist =true;
        return position; // return the position
    }
    else if (mark == 1){
        exist =false;
        return position + 1; // insert after the position
    }
    else{  // mark == -1.
        exist = false;
        return position;  // insert before the position.
    }
}

// --------------------- Use merge sort algorithm to sort sys_quotient_state ------------------------------
sys_quotient_state combine_sys_quotient_state(const sys_quotient_state & s1, const sys_quotient_state & s2){
    int i;
    int size2= s2.sysxindex.size();
    sys_quotient_state s3=s1;
    for(i=0;i<size2;i++){
        s3.sysxindex.push_back(s2.sysxindex[i]);
        s3.xindex.push_back(s2.xindex[i]);
    }
    return s3;
}


vector<sys_quotient_state> merge_sort(vector<sys_quotient_state> & list1, vector<sys_quotient_state> & list2){
    int j;
    int size1= list1.size();
    int size2= list2.size();
    int list1_index=0;
    int list2_index=0;
    int mark;
    if (size1==0) return list2;
    if(size2==0) return list1;
    vector<sys_quotient_state> list_after_merge;
    sys_quotient_state s3=list1[0];

    while(list1_index< size1 and list2_index<size2){
        mark= compare_vmode(list1[list1_index].vmode1,list1[list1_index].vmode2,list2[list2_index].vmode1,list2[list2_index].vmode2);
        if(mark>0){
            // list2 's element should in front of list1's element.
            list_after_merge.push_back(list2[list2_index]);
            list2_index++;
        }
        else if (mark<0){
            list_after_merge.push_back(list1[list1_index]);
            list1_index++;
        }
        else if (mark==0){
            // combine two sys_quotient_state
            s3= combine_sys_quotient_state(list1[list1_index],list2[list2_index]);
            list_after_merge.push_back(s3);
            list1_index++;
            list2_index++;
        }
    }
    int i;
    // push_back element left in list
    if(list1_index < size1){
        for(i=list1_index;i<size1;i++){
            list_after_merge.push_back(list1[i]);
        }
    }
    if(list2_index< size2){
        for(i=list2_index;i<size2;i++){
            list_after_merge.push_back(list2[i]);
        }
    }
    return list_after_merge;
}

vector<sys_quotient_state> merge_sort_list(vector<sys_quotient_state> & list){
    int list_size= list.size();
    int i;
    vector<vector<sys_quotient_state>> * new_list_ptr;
    vector<vector<sys_quotient_state>> * old_list_ptr;
    vector<vector<sys_quotient_state>> *list_ptr_3;
    vector<vector<sys_quotient_state>> old_list;
    for(i=0;i<list_size;i++){
        vector<sys_quotient_state> list_element;
        list_element.push_back(list[i]);
        old_list.push_back(list_element);
    }
    old_list_ptr= & (old_list);
    vector<vector<sys_quotient_state>> new_list;
    new_list_ptr = & new_list;

    while (list_size>1){
        vector<sys_quotient_state> merge_result;
        for(i=0;i+1<list_size;i=i+2){
            merge_result= merge_sort((*old_list_ptr)[i],(*old_list_ptr)[i+1]);
            (*new_list_ptr).push_back(merge_result);
        }
        if(list_size % 2 ==1){
            (*new_list_ptr).push_back((*old_list_ptr)[list_size-1]);  // add last list if size of list is odd.
        }
        (*old_list_ptr).clear();
        (*old_list_ptr).shrink_to_fit();   // free the space
        // exchange list
        list_ptr_3 = &(*old_list_ptr);
        old_list_ptr =  &(*new_list_ptr);  // use new list.
        new_list_ptr= & (*list_ptr_3);

        list_size= (list_size +1)/2;  // update size of list after merge sort.
    }
    vector <sys_quotient_state> final_list;
    if(list_size!=0){
        final_list = (*old_list_ptr)[0];
    }

    return final_list;
}

// ---------- end of merge_sort algorithm ----------------------------------------------

// save sys quotient state information.
void save_sys_quotient_state_data(const vector<sys_quotient_state> & list, string path){
    //  format to store data :  list_size  moddim
    // every line is one sys_quotient_state struct:
    // vmode1, vmode2, xindex_length (sysindex_length), xindex,  sysxindex
    ofstream save(path + "sys_quotient_state.txt");
    int size= list.size();
    int moddim = list[0].vmode1.size();
    int xlistsize;
    int i,j;
    save<< size << " "<< moddim<<endl;
    for(i=0;i<size;i++){
        // output vmode1:
        for(j=0;j<moddim;j++){
            save<< list[i].vmode1[j]<< " ";
        }
        for(j=0;j<moddim;j++){
            save<< list[i].vmode2[j]<<" ";
        }
        save<< list[i].dindex1<<" "<<list[i].dindex2<<" ";
        xlistsize= list[i].xindex.size();
        save<< xlistsize << " ";
        for(j=0;j<xlistsize;j++){
            save<< list[i].xindex[j]<<" ";
        }
        for(j=0;j<xlistsize;j++){
            save<<list[i].sysxindex[j]<<" ";
        }
        save<<endl;
    }
    save.close();
}

void load_sys_quotient_state_data( vector <sys_quotient_state> &list, string path){
    //  format to load data :  list_size  moddim
    // every line is one sys_quotient_state struct:
    // vmode1, vmode2, xindex_length (sysindex_length), xindex,  sysxindex
    ifstream load (path +  "sys_quotient_state.txt");  // load data from sys_quotient_state.txt file.
    int size;
    int moddim;
    int xlistsize;
    int i,j;
    int x;
    int dindex1, dindex2;
    string ss;
    load>> size >> moddim;
    std::getline(load, ss);
    for(i=0;i<size;i++){
        vector <int> vmode1;
        vector <int> vmode2;
        for(j=0;j<moddim;j++){
            load >> x;
            vmode1.push_back(x);
        }
        for(j=0;j<moddim;j++){
            load>> x;
            vmode2.push_back(x);
        }
        load>> dindex1 >> dindex2;
        sys_quotient_state sq (vmode1,vmode2,dindex1,dindex2);
        load>> xlistsize;
        for(j=0;j<xlistsize;j++){
            load >> x;
            sq.xindex.push_back (x);
        }
        for(j=0;j<xlistsize;j++){
            load>> x;
            sq.sysxindex.push_back(x);
        }
        std::getline(load, ss);
        list.push_back(sq);
    }
    load.close();
}