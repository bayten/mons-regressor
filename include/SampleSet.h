/* Copyright 2017 Baytekov Nikita */
#ifndef INCLUDE_SAMPLESET_H_
#define INCLUDE_SAMPLESET_H_

#include <algorithm>
#include <vector>
#include "default_types.h"


template<typename S, typename T>  // S - type of features, T - type of answers
class GroupSamples {
 private:
    Mat<S> objs;
    T group_tag;

 public:
    GroupSamples() : objs(0), group_tag() {}
    GroupSamples(Mat<S> init_objs, T init_tag) :
                objs(init_objs), group_tag(init_tag) {}
    GroupSamples(Vec<S> init_vec, T init_tag) :
                objs(Mat<S>(init_vec)), group_tag(init_tag) {}
    ~GroupSamples() {}

    void shuffle() { std::random_shuffle(&(objs[0]), &(objs[-1])+1); }
    void slice_rand(int slice_num, Mat<S>* slice_x, Vec<T>* slice_y) const;
    void append(const Vec<S>& X) { objs.hadd(X); }
    void append(const Mat<S>& X);

    int get_size() const { return objs.get_sx(); }
    T get_tag() const { return group_tag; }
    const Mat<S>& get_objs() const { return objs; }

    Vec<S>& operator[](int index) { return objs[index]; }
    const Vec<S>& operator[](int index) const { return objs[index]; }

    template<typename U, typename V>
    friend std::ostream& operator<<(std::ostream& os, const GroupSamples<U, V>& gr_samps);
};

template<typename S, typename T>  // S - type of features, T - type of answers
class SampleSet {
 private:
    Vec< GroupSamples<S, T> > groups;

 public:
    explicit SampleSet(int size = 0) : groups(size) {}  // default constructor
    explicit SampleSet(Vec< GroupSamples<S, T> > init_groups);
    SampleSet(const SampleSet<S, T>& copy_obj) :  // copy constructor
                groups(copy_obj.groups) {}

    explicit SampleSet(SampleSet<S, T>&& move_obj) :  // move constructor
                groups(std::move(move_obj.get_groups())) {}
    ~SampleSet() {}

    void append(const Mat<S>& X, const Vec<T>& y);
    void shuffle();
    int get_group_num() const { return groups.get_size(); }
    int get_total_size() const;
    Vec< GroupSamples<S, T> >& get_groups() { return groups; }
    const Vec< GroupSamples<S, T> >& get_groups() const { return groups; }
    const GroupSamples<S, T>& get_group(T index_tag) const;
    SampleSet<S, T> get_antigroup(T index_tag) const;
    bool delete_group(T index_tag);

    SampleSet<S, T>& operator= (const SampleSet<S, T>& copy_obj);  // copy assignment
    Vec<S>& operator[](int abs_index);
    const Vec<S>& operator[](int abs_index) const;

    template<typename U, typename V>
    friend std::ostream& operator<<(std::ostream& os, const SampleSet<U, V>& sset);
};


template<typename S, typename T>
void GroupSamples<S, T>::slice_rand(int slice_num, Mat<S>* slice_x, Vec<T>* slice_y) const {
    std::vector<int>idx(objs.get_sx());
    std::iota(std::begin(idx), std::end(idx), 0);
    std::random_shuffle(std::begin(idx), std::end(idx));

    *slice_x = Mat<S>(slice_num, objs.get_sy());
    *slice_y = Vec<T>(slice_num);

    for (int i = 0; i < slice_num; i++) {
        (*slice_x)[i] = objs[i];
        (*slice_y)[i] = group_tag;
    }
}

template<typename S, typename T>
void GroupSamples<S, T>::append(const Mat<S>& X) {
    int mat_sx = X.get_sx();

    for (int i = 0; i < mat_sx; i++)
        objs.append(X[i]);
}

template<typename U, typename V>
std::ostream& operator<<(std::ostream& os, const GroupSamples<U, V>& gr_samps) {
    int gsamps_size = gr_samps.objs.get_sx();
    os << "GSamps{" << gsamps_size << ", tag:" << gr_samps.group_tag << "}:[" << gr_samps.objs;
    return os;
}

template<typename S, typename T>
SampleSet<S, T>::SampleSet(Vec< GroupSamples<S, T> > init_groups) :
        groups(init_groups) {}

template<typename S, typename T>
void SampleSet<S, T>::append(const Mat<S>& X, const Vec<T>& y) {
    int obj_num = y.get_size();
    int i = 0;

    if (!groups.get_size()) {  // if appending for the first time...
        // LOG_(trace) << "Appending for the first time...";
        groups.append(GroupSamples<S, T>(X[0], y[0]));
        i++;
    }

    for (; i < obj_num; i++) {
        bool was_found = 0;

        for (int j = 0; j < groups.get_size(); j++) {  // searching within existing groupes
            if (groups[j].get_tag() == y[i]) {
                groups[j].append(X[i]);
                was_found = 1;
                break;
            }
        }
        if (!was_found)  // otherwise adding new group container =)
            groups.append(GroupSamples<S, T>(X[i], y[i]));
    }
    // LOG_(trace) << "Matrix of objects was appended: " << (*this);
}

template<typename S, typename T>
void SampleSet<S, T>::shuffle() {
    int groups_size = groups.get_size();
    for (int i = 0; i < groups_size; i++)
        groups[i].shuffle();
}

template<typename S, typename T>
int SampleSet<S, T>::get_total_size() const {
    int groups_num = groups.get_size();
    int total_size = 0;
    for (int i = 0; i < groups_num; i++)
        total_size += groups[i].get_size();

    return total_size;
}

template<typename S, typename T>
const GroupSamples<S, T>& SampleSet<S, T>::get_group(T index_tag) const {
    int groups_size = groups.get_size();
    for (int i = 0; i < groups_size; i++)
        if (groups[i].get_tag() == index_tag)
            return groups[i];

    return groups[0];
}

template<typename S, typename T>
SampleSet<S, T> SampleSet<S, T>::get_antigroup(T index_tag) const {
    SampleSet<S, T> new_sample_set(*this);
    new_sample_set.delete_group(index_tag);
    return new_sample_set;
}

template<typename S, typename T>
bool SampleSet<S, T>::delete_group(T index_tag) {
    int groups_size = groups.get_size();
    for (int i = 0; i < groups_size; i++) {
        if (groups[i].get_tag() == index_tag) {
            groups.erase(i);
            return 1;
        }
    }
    LOG_(warning) << "No group samples with tag " << index_tag << " were found.";
    return 0;
}

template<typename S, typename T>
Vec<S>& SampleSet<S, T>::operator[](int abs_index) {
    int groups_num = groups.get_size();
    for (int i = 0; i < groups_num; i++) {
        if (abs_index < groups[i].get_size())
            return groups[i][abs_index];

        abs_index -= groups[i].get_size();
    }

    return groups[-1][-1];
}

template<typename S, typename T>
const Vec<S>& SampleSet<S, T>::operator[](int abs_index) const {
    int groups_num = groups.get_size();
    for (int i = 0; i < groups_num; i++) {
        if (abs_index < groups[i].get_size())
            return groups[i][abs_index];

        abs_index -= groups[i].get_size();
    }

    return groups[-1][-1];
}

template<typename S, typename T>
SampleSet<S, T>& SampleSet<S, T>::operator= (const SampleSet<S, T>& copy_obj) {
    groups = copy_obj.get_groups();
    return (*this);
}

template<typename U, typename V>
std::ostream& operator<<(std::ostream& os, const SampleSet<U, V>& sset) {
    std::stringstream buffer;
    int sset_size = sset.groups.get_size();
    int total_size = sset.get_total_size();
    buffer << "SampleSet{" << total_size << "; " << sset_size << "}:[";

    if (sset_size < 1) {
        buffer << "nullptr";
    } else if (sset_size < 10) {
        buffer << sset.groups[0];
        for (int i = 1; i < sset_size; i++)
            buffer << ", " << sset.groups[i] << std::endl;
    } else {
        buffer << sset.groups[0] << "," << std::endl;
        buffer << sset.groups[1] << "," << std::endl;
        buffer << sset.groups[2] << "," << std::endl << "..., " << std::endl;
        buffer << sset.groups[sset_size-3] << "," << std::endl;
        buffer << sset.groups[sset_size-2] << "," << std::endl;
        buffer << sset.groups[sset_size-1] << "," << std::endl;
    }
    buffer << "]";
    os << buffer.str();
    return os;
}

#endif  // INCLUDE_SAMPLESET_H_
