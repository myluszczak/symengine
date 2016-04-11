/**
 *  \file series_generic.h
 *  Class for univariate series.
 *
 **/
#ifndef SYMENGINE_SERIES_GENERIC_H
#define SYMENGINE_SERIES_GENERIC_H

#include <list>
#include <symengine/dict.h>
#include <symengine/polynomial.h>
#include <symengine/series.h>
//#include <symengine/series_visitor.h>

namespace SymEngine
{
//! UnivariateSeries Class
class UnivariateSeries : public Basic
{
public:
    //! `var_` : Variable of the UnivariateSeries
    //! `poly_` : holds the UnivariateIntPolynomial of the series
    //! `prec_` : precision of the UnivariateSeries, i.e. self = poly +
    //! O(var^prec)
    // UnivariateSeries 1 + 2*x + x**2 + O(x**5) has dict_ = {{0, 1}, {1, 2},
    // {2, 1}} with var_ = "x" and prec_ = 5
    RCP<const Symbol> var_;
    RCP<const UnivariateIntPolynomial> poly_;
    unsigned int prec_;

public:
    IMPLEMENT_TYPEID(UNIVARIATESERIES)
    //! Constructor of UnivariateSeries class
    UnivariateSeries(const RCP<const Symbol> &var,
                     const unsigned int &precision,
                     const RCP<const UnivariateIntPolynomial> &poly);
    UnivariateSeries(const RCP<const Symbol> &var,
                     const unsigned int &precision, const unsigned int &max_exp,
                     map_uint_mpz &&dict);
    UnivariateSeries(const RCP<const Symbol> &var,
                     const unsigned int &precision, const map_uint_mpz &dict);
    //! Constructor using a dense vector of integer_class coefficients
    UnivariateSeries(const RCP<const Symbol> &var,
                     const unsigned int &precision,
                     const std::vector<integer_class> &v);

    static RCP<const UnivariateSeries>
    create(const RCP<const Symbol> &var, const unsigned int &prec,
           const std::vector<integer_class> &v)
    {
        return make_rcp<const UnivariateSeries>(var, prec, v);
    }

    //! \return true if canonical
    bool is_canonical(const UnivariateIntPolynomial &,
                      const unsigned int &) const;
    //! \return size of the hash
    std::size_t __hash__() const;
    /*! Equality comparator
     * \param o - Object to be compared with
     * \return whether the 2 objects are equal
     * */

    bool operator==(const UnivariateSeries &o) const
    {
        return (eq(*var_, *o.var_) and poly_->__eq__(*o.poly_)
                and prec_ == o.prec_);
    }
    bool __eq__(const Basic &o) const;
    int compare(const Basic &o) const;

    std::string __str__() const;
    virtual vec_basic get_args() const
    {
        return {};
    }
};

inline RCP<const UnivariateSeries> univariate_series(RCP<const Symbol> i,
                                                     unsigned int prec,
                                                     const map_uint_mpz &dict)
{
    return make_rcp<const UnivariateSeries>(i, prec, dict);
}

RCP<const UnivariateSeries> add_uni_series(const UnivariateSeries &a,
                                           const UnivariateSeries &b);
RCP<const UnivariateSeries> neg_uni_series(const UnivariateSeries &a);
RCP<const UnivariateSeries> sub_uni_series(const UnivariateSeries &a,
                                           const UnivariateSeries &b);
RCP<const UnivariateSeries> mul_uni_series(const UnivariateSeries &a,
                                           const UnivariateSeries &b);

class MultivariateSeries : public SeriesBase<MultivariateExprPolynomial, Expression, MultivariateSeries>{
public:
    set_sym vars_;
    umap_sym_uint precs_;
public:
    IMPLEMENT_TYPEID(MULTIVARIATESERIES);
    MultivariateSeries(const MultivariateExprPolynomial p, const set_sym vars, const umap_sym_uint precs_) :
        SeriesBase(std::move(p),"",0), vars_{std::move(vars)}, precs_{std::move(precs_)}
    {
    }

    inline bool __eq__(const Basic &o) const
    {
        return (is_a<MultivariateSeries>(o) and vars_ == static_cast<const MultivariateSeries &>(o).vars_
                and p_ == static_cast<const MultivariateSeries &>(o).p_
                and precs_ == static_cast<const MultivariateSeries &>(o).precs_);
    }
    
    std::size_t __hash__() const {
        //DUMMY
        return 0;
    }
    
    int compare(const Basic &o) const {
        //DUMMY
        return 0;
    }
   
    RCP<const Basic> as_basic() const {
        //DUMMY
        RCP<const Basic> b; 
        return b; 
    }

    umap_int_basic as_dict() const {
        //DUMMY
        umap_int_basic m;
        return m;
    }

    RCP<const Basic> get_coeff(int) const {
        //DUMMY
        return make_rcp<const Integer>(0);
    }

    RCP<const Number> add(const Number &other) const
    {/*
        if (is_a<MultivariateSeries>(other)) {
            const MultivariateSeries &o = static_cast<const MultivariateSeries &>(other);
            umap_sym_uint precs;
            //take the minimum precision for each variable
            for (auto bucket : precs_) {
                precs.insert(std::pair<RCP<const Symbol>, unsigned int>(bucket.first,bucket.second) );
            }
            for (auto bucket : o.precs_) {
                if (precs.find(bucket.first) != precs.end()) {
                    precs.insert(std::pair<RCP<const Symbol>, unsigned int>(bucket.first,bucket.second));
                } else {
                    if (precs.find(bucket.first)->second > bucket.second) 
                        precs.find(bucket.first)->second = bucket.second;
                }
            }
            return make_rcp<MultivariateSeries>(MultivariateExprPolynomial(p_ + o.p_), vars_, precs);
        } else if (other.get_type_code() < MultivariateSeries::type_code_id) {
            MultivariateExprPolynomial p = MultivariateSeries::series(other.rcp_from_this(), vars_, precs_)->p_;
            return make_rcp<MultivariateSeries>(MultivariateExprPolynomial(p_ + p), vars_, precs_);
        } else {
            return other.add(*this);
        }*/
        RCP<const Number> n;
        return n;
    }

    RCP<const Number> mul(const Number &other) const
    {/*
        if (is_a<MultivariateSeries>(other)) {
            const MultivariateSeries &o = static_cast<const MultivariateSeries &>(other);
            umap_sym_uint precs;
            for (auto bucket : precs_) {
                precs.insert(std::pair<RCP<const Symbol>, unsigned int>(bucket.first,bucket.second) );
            }
            for (auto bucket : o.precs_) {
                if (precs.find(bucket.first) != precs.end()) {
                    precs.insert(std::pair<RCP<const Symbol>, unsigned int>(bucket.first,bucket.second));
                } else {
                    if (precs.find(bucket.first)->second > bucket.second) 
                        precs.find(bucket.first)->second = bucket.second;
                }
            }
            return make_rcp<MultivariateSeries>(MultivariateSeries::mul(p_, o.p_, precs), vars_, precs);
        } else if (other.get_type_code() < MultivariateSeries::type_code_id) {
            MultivariateExprPolynomial p = MultivariateSeries::series(other.rcp_from_this(), vars_, precs)->p_;
            return make_rcp<MultivariateSeries>(MultivariateSeries::mul(p_, p, precs), vars_, precs);
        } else {
            return other.mul(*this);
        }*/
        RCP<const Number> n;
        return n;
    }

    RCP<const Number> pow(const Number &other) const
    {/*
        umap_sym_uint precs;
        for (auto bucket : precs_) {
            precs.insert(std::pair<RCP<const Symbol>, unsigned int>(bucket.first,bucket.second) );
        }
        MultivariateExprPolynomial p;
        if (is_a<MultivariateSeries>(other)) {
            const MultivariateSeries &o = static_cast<const MultivariateSeries &>(other);

            for (auto bucket : o.precs_) {
                if (precs.find(bucket.first) != precs.end()) {
                    precs.insert(std::pair<RCP<const Symbol>, unsigned int>(bucket.first,bucket.second));
                } else {
                    if (precs.find(bucket.first)->second > bucket.second) 
                        precs.find(bucket.first)->second = bucket.second;
                }
            }
            p = o.p_;
        } else if (is_a<Integer>(other)) {
            if (other.is_negative()) {
                p = MultivariateSeries::pow(
                    p_, (static_cast<const Integer &>(other).neg()->as_int()),
                    precs);
                p = MultivariateSeries::series_invert(p, MultivariateSeries::vars(vars_), precs);
                return make_rcp<MultivariateSeries>(p, vars_, precs);
            }
            p = MultivariateSeries::pow(p_, (static_cast<const Integer &>(other).as_int()),
                            precs);
            return make_rcp<MultivariateSeries>(p, vars_, precs);
        } else if (other.get_type_code() < Series::type_code_id) {
            p = MultivariateSeries::series(other.rcp_from_this(), vars_, precs)->p_;
        } else {
            return other.rpow(*this);
        }
        p = MultivariateSeries::series_exp(
            MultivariateExprPolynomial(p * MultivariateSeries::series_log(p_, MultivariateSeries::vars(vars_), precs)),
            MultivariateSeries::vars(vars_), precs);
        return make_rcp<MultivariateSeries>(p, vars_, precs);*/
        RCP<const Number> n;
        return n;
    }

    RCP<const Number> rpow(const Number &other) const
    {/*
        if (other.get_type_code() < MultivariateSeries::type_code_id) {
            MultivariatePolynomial p = MultivariateSeries::series(other.rcp_from_this(), vars_, precs)->p_;
            p = MultivariateSeries::series_exp(
                MultivariateExprPolynomial(p_ * MultivariateSeries::series_log(p, MultivariateSeries::vars(vars_), precs)),
                MultivariateSeries::vars(vars_), precs);
            return make_rcp<MultivariateSeries>(p, vars_, prec);
        } else {
            throw std::runtime_error("Unknown type");
        }*/
        RCP<const Number> n;
        return n;
    }

    RCP <const MultivariateSeries> series(RCP<const Basic> b, set_sym vars, umap_sym_uint precs) const {
        //SeriesVisitor<MultivariateExprPolynomial, Expression, MultivariateSeries>(MultivariatePolynomial({},{},{}), "", 0);
        //return visitor.series(b);
        RCP<const MultivariateSeries> s;
        return s;
    }

 
};
 
} // SymEngine
#endif
