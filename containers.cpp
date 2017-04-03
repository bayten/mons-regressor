#include <algorithm>


template<class T>  // type of features
struct ElClass_t 
{
  public:
    Vec_t<int> cols;
    Vec_t<T> vals;

    ElClass_t(int rank) : cols(rank), vals(rank)
    {
        for(int i = 0; i < rank; i++)
            cols[i] = -1;
    }

    ElClass_t(int rank, int init_cols[], T init_vals[]) : cols(rank), vals(rank)
    {
        for(int i = 0; i < rank; i++){
            cols[i] = init_cols[i];
            vals[i] = init_vals[i];
        }
    }

    ~ElClass_t(){};
};


template<class T>  
struct ElColl_t
{
    Vec_t< ElClass_t<T> > ecs;

    ElColl_t(int size) : ecs(size){};

    ElColl_t(int size, ElClass_t<T> init_ecs[]) : ecs(size)
    {
        for(int i = 0; i < size; i++)
            ecs[i] = init_ecs[i];
    }

    ~ElColl_t(){};
};


template<class T1, class T2> 
class ClassSamples_t
{
  public:
    Mat_t<T1> objs;
    T2 class_tag;

    ClassSamples_t(Mat_t<T1> init_objs, T2 init_tag) : objs(init_objs), class_tag(init_tag){};
    ~ClassSamples_t(){};

    void tearup(int tear_num, Mat_t<T1> & tornX, Vec_t<T2> & torny)
    {
        std::random_shuffle( &(objs[0]), &(objs[-1]) );
        tornX = objs.cutslice(0, tear_num);
        torny = Vec_t<T2>(tear_num);
        for(int j = 0; j < tear_num; j++)
            torny[j] = class_tag;
    }
};


template<class T1, class T2>
class SampleSet_t
{
  public:
    Vec_t< ClassSamples_t<T1, T2> > samples;

    SampleSet_t(int size=0) : samples(size) {};
    ~SampleSet_t(){};

    void append(Mat_t<T1> X, Vec_t<T2> y)
    {
        int obj_num = y.sz;
        int i = 0;
        bool was_found = 0;

        if(!samples.sz){  // if appending for the first time...
            samples.append(ClassSamples_t<T1, T2>(X[0], y[0]));
            i++;
        }

        for(; i < obj_num; i++){
            was_found = 0;

            for(int j = 0; j < samples.sz; j++){  // searching within existing classes
                if(samples[j].class_tag == y[i]){
                    samples[j].objs.append(X[i]);
                    was_found = 1;
                    break;
                }
            }
            if(!was_found){  // otherwise adding new class container =)
                samples.append(ClassSamples_t<T1, T2>(X[i], y[i]));
            }
        }
    }

    void shuffle()
    {
        for(int i = 0; i < samples.sz; i++){
            std::random_shuffle(&(samples[i].obj[0]), &(samples[i].obj[-1]));
        }
    }


    ClassSamples_t<T1, T2> & operator[](T2 index_tag)
    {
        for(int i = 0; i < samples.sz; i++)
            if(samples[i].class_tag == index_tag)
                return samples[i];

        return sample[0];
    }
};
