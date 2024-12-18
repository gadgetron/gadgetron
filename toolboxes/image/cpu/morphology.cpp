/**
\file          morphology.cpp
\brief      Implement morphological image processing related functions.
\author     Hui Xue
*/

#include "morphology.h"
#include <iostream>
#include <stack>
#include <cmath>

namespace Gadgetron
{

template <typename T>
void region_growing_2d(const hoNDArray<T>& input, T object_value, hoNDArray<unsigned int>& label_array, size_t x, size_t y, unsigned int label, bool is_8_connected)
{
    try
    {
        typedef std::pair<long long, long long> PtType;

        std::stack< std::pair<long long, long long> > ss;

        PtType seed;
        seed.first = x;
        seed.second = y;

        ss.push(seed);

        std::vector<long long> dx, dy;
        if (is_8_connected)
        {
            dx.resize(8);
            dy.resize(8);

            dx[0] = -1;
            dy[0] = -1;

            dx[1] = 0;
            dy[1] = -1;

            dx[2] = 1;
            dy[2] = -1;

            dx[3] = -1;
            dy[3] = 0;

            dx[4] = 1;
            dy[4] = 0;

            dx[5] = -1;
            dy[5] = 1;

            dx[6] = 0;
            dy[6] = 1;

            dx[7] = 1;
            dy[7] = 1;
        }
        else
        {
            dx.resize(4);
            dy.resize(4);

            dx[0] = 0;
            dy[0] = -1;

            dx[1] = -1;
            dy[1] = 0;

            dx[2] = 1;
            dy[2] = 0;

            dx[3] = 0;
            dy[3] = 1;
        }

        size_t n;
        size_t numNeighbor = dx.size();

        while (!ss.empty())
        {
            PtType v = ss.top();
            ss.pop();

            label_array(v.first, v.second) = label;
            // GDEBUG_STREAM("label point : " << v.first << " " << v.second);

            long long cx, cy;

            for (n = 0; n < numNeighbor; n++)
            {
                cx = v.first + dx[n];
                cy = v.second + dy[n];

                if (input.point_in_range(cx, cy))
                {
                    if (label_array(cx, cy) == 0 && std::abs(input(cx, cy) - object_value) < FLT_EPSILON)
                    {
                        ss.push(PtType(cx, cy));
                    }
                }
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors happened in region_growing_2d(...) ... ");
    }
}

template void region_growing_2d(const hoNDArray<int>& input, int object_value, hoNDArray<unsigned int>& label_array, size_t x, size_t y, unsigned int label, bool is_8_connected);
template void region_growing_2d(const hoNDArray<float>& input, float object_value, hoNDArray<unsigned int>& label_array, size_t x, size_t y, unsigned int label, bool is_8_connected);
template void region_growing_2d(const hoNDArray<double>& input, double object_value, hoNDArray<unsigned int>& label_array, size_t x, size_t y, unsigned int label, bool is_8_connected);

// --------------------------------------------------------------------------------------------

template <typename T>
void bwlabel_2d(const hoNDArray<T>& input, T object_value, hoNDArray<unsigned int>& label, bool is_8_connected)
{
    try
    {
        size_t COL = input.get_size(0);
        size_t ROW = input.get_size(1);
        size_t num = input.get_number_of_elements();

        label.create(COL, ROW);
        Gadgetron::clear(label);

        size_t currLabel = 1;

        size_t r, c;
        for (r = 1; r < ROW - 1; r++)
        {
            for (c = 1; c < COL - 1; c++)
            {
                if (std::abs(input(c, r) - object_value) < FLT_EPSILON)
                {
                    if (label(c, r) == 0)
                    {
                        Gadgetron::region_growing_2d(input, object_value, label, c, r, currLabel, is_8_connected);
                        currLabel++;
                    }
                }
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors happened in bwlabel_2d(...) ... ");
    }
}

template void bwlabel_2d(const hoNDArray<int>& input, int object_value, hoNDArray<unsigned int>& label, bool is_8_connected);
template void bwlabel_2d(const hoNDArray<float>& input, float object_value, hoNDArray<unsigned int>& label, bool is_8_connected);
template void bwlabel_2d(const hoNDArray<double>& input, double object_value, hoNDArray<unsigned int>& label, bool is_8_connected);

// --------------------------------------------------------------------------------------------

template <typename T>
void bwlabel_clean_fore_and_background(const hoNDArray<T>& input, T object_value, T bg_value, size_t obj_thres, size_t bg_thres, bool is_8_connected, hoNDArray<T>& output)
{
    try
    {
        output = input;

        size_t RO = input.get_size(0);
        size_t E1 = input.get_size(1);

        // first, clean background
        hoNDArray<unsigned int> label;
        label.clear();

        Gadgetron::bwlabel_2d(input, bg_value, label, is_8_connected);

        std::vector<unsigned int> labels, areas;
        Gadgetron::bwlabel_area_2d(label, labels, areas);

        size_t nn, ro, e1;
        size_t num_areas = areas.size();
        for (nn = 0; nn<num_areas; nn++)
        {
            if (areas[nn]<bg_thres)
            {
                for (e1 = 0; e1 < E1; e1++)
                {
                    for (ro = 0; ro < RO; ro++)
                    {
                        if (label(ro, e1) == labels[nn])
                        {
                            output(ro, e1) = object_value;
                        }
                    }
                }
            }
        }

        // clean forground
        label.clear();
        Gadgetron::bwlabel_2d(output, object_value, label, is_8_connected);
        Gadgetron::bwlabel_area_2d(label, labels, areas);

        num_areas = areas.size();
        for (nn = 0; nn<num_areas; nn++)
        {
            if (areas[nn]<obj_thres)
            {
                for (e1 = 0; e1 < E1; e1++)
                {
                    for (ro = 0; ro < RO; ro++)
                    {
                        if (label(ro, e1) == labels[nn])
                        {
                            output(ro, e1) = bg_value;
                        }
                    }
                }
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors happened in bwlabel_clean_fore_and_background(...) ... ");
    }
}

template void bwlabel_clean_fore_and_background(const hoNDArray<int>& input, int object_value, int bg_value, size_t obj_thres, size_t bg_size, bool is_8_connected, hoNDArray<int>& output);
template void bwlabel_clean_fore_and_background(const hoNDArray<float>& input, float object_value, float bg_value, size_t obj_thres, size_t bg_size, bool is_8_connected, hoNDArray<float>& output);
template void bwlabel_clean_fore_and_background(const hoNDArray<double>& input, double object_value, double bg_value, size_t obj_thres, size_t bg_size, bool is_8_connected, hoNDArray<double>& output);

// --------------------------------------------------------------------------------------------

void bwlabel_area_2d(const hoNDArray<unsigned int>& label_array, std::vector<unsigned int>& labels, std::vector<unsigned int>& areas)
{
    try
    {
        labels.clear();
        areas.clear();

        size_t num = label_array.get_number_of_elements();

        size_t n;
        for (n = 0; n < num; n++)
        {
            unsigned int v = label_array(n);
            if (v>0)
            {
                size_t numLabel = labels.size();

                bool labelExist = false;
                for (size_t ll = 0; ll < numLabel; ll++)
                {
                    if (labels[ll] == v)
                    {
                        labelExist = true;
                        break;
                    }
                }

                if (!labelExist)
                {
                    labels.push_back(v);
                }
            }
        }

        size_t numLabels = labels.size();

        if (numLabels == 0) return;

        areas.resize(numLabels, 0);

        for (n = 0; n < num; n++)
        {
            unsigned int v = label_array(n);
            if (v>0)
            {
                size_t l;
                for (l = 0; l < numLabels; l++)
                {
                    if (v == labels[l]) break;
                }

                if (l<numLabels)
                {
                    areas[l]++;
                }
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors happened in bwlabel_area_2d(...) ... ");
    }
}

}
