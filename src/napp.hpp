/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 */

#ifndef _NAMED_ARGUMENTS_NA_HPP
#define _NAMED_ARGUMENTS_NA_HPP 1

#include <tuple>
#include <any>
#include <memory>

#if 0
#if 0
template<typename U, typename... T>
constexpr bool contains(std::tuple<T...>) {
    return (std::is_same_v<U, T> || ...);
}

#else
template <typename T, typename Tuple>
struct has_type;

template <typename T, typename... Us>
struct has_type<T, std::tuple<Us...>> : std::disjunction<std::is_same<T, Us>...> {};
#endif
#endif


namespace NA {

namespace detail
{

struct owner_object {};
struct reference_object {};

template <typename T, typename tag,typename ObjectType>
struct type;

template <typename T>
struct ownership_object : std::conditional<
    std::is_lvalue_reference_v<T>,
    reference_object,
    owner_object>
{};

template <typename T,typename tag>
struct infer_type
{
    using object_type = typename ownership_object<T>::type;
    using value_type = std::remove_reference_t<T>;
    using type = NA::detail::type<value_type,tag,object_type>;
};

template <typename T,typename tag>
using infer_type_t = typename infer_type<T,tag>::type;


struct named_argument_base {};





template <typename T, typename tag>
struct type<T,tag,owner_object> : named_argument_base
{
    static_assert(!std::is_reference_v<T>,
                  "T must not be a reference. Rather set the category!");
    using tag_type = tag;
    using value_type = T;

    template <typename T2>
    using apply_t = type<T,tag,owner_object>;

    template <typename T2>
    using required_as_t = typename detail::infer_type_t<T2,tag>;


    type(type &&) = default;
    type(type const&) = delete;

    template <typename V, std::enable_if_t< std::is_constructible_v<T,V&&> , bool> = true>
    constexpr type(V && t) : M_value(std::forward<V>(t)) {}

    // build from another owner_object
    template <typename T2, std::enable_if_t< std::is_constructible_v<T,T2&&> , bool> = true>
    constexpr type( type<T2,tag,owner_object> && t ) : M_value(std::move(t).value()) {}


    template <typename T2>
    constexpr type( type<T2,tag,reference_object> && t ) : M_value(std::move(t).value()) { /*std::cout << "owner_object-2-" << std::endl;*/ } // TODO here no apply copy
    //type( type<T2,tag,reference_object> && t ) : M_value(std::forward<type<T2,tag,reference_object>>(t).value()) { std::cout << "owner_object-2-" << std::endl; } // TODO here no apply copy

    constexpr T & value() & { return M_value; }
    constexpr T const & value() const & { return M_value; }
    constexpr T && value() && { return std::move( M_value ); }

private :
    T M_value;
};

template <typename T, typename tag>
struct type<T,tag,reference_object> : named_argument_base
{
    using tag_type = tag;
    using value_type = T;

    template <typename T2>
    using apply_t = type<T,tag,reference_object>;

    template <typename T2>
    using required_as_t = typename detail::infer_type_t<T2,tag>;

    constexpr type(T & t) : M_value(t) {}

    type(type &&) = default;//delete;
    type(type const&) = delete;

    template <typename V,std::enable_if_t< std::is_lvalue_reference_v<V&&> && std::is_base_of_v<std::decay_t<T>,std::decay_t<V>> , bool> = true>
    constexpr type(V && t) : M_value(std::forward<V>(t)) {}

    template <typename V,std::enable_if_t< std::is_rvalue_reference_v<V&&>   && !std::is_base_of_v<named_argument_base,std::decay_t<V> >      /*&& !std::is_same_v<std::decay_t<T>,std::decay_t<V>>*/, bool> = true>
    constexpr type(V && t) : M_tmp(std::make_shared<T>(std::forward<V>(t))), M_value( *std::any_cast<std::shared_ptr<T>&>(M_tmp).get() ) {}

    template <typename V,std::enable_if_t< std::is_lvalue_reference_v<V&&> && !std::is_base_of_v<std::decay_t<T>,std::decay_t<V>>, bool> = true>
    constexpr type(V && t) : M_tmp(std::make_shared<T>(std::forward<V>(t))), M_value( *std::any_cast<std::shared_ptr<T>&>(M_tmp).get() ) {}


    template <typename T2, typename tag2,typename ObjectType2> friend struct type;


    template <typename T2>
    constexpr type( type<T2,tag,owner_object> && t ) :
        M_tmp(std::make_shared<T>(std::forward<type<T2,tag,owner_object>>(t).value())),
        M_value( *std::any_cast<std::shared_ptr<T>&>(M_tmp).get() ) {}


    // TODO HERE : maybe some case where we can avoid the copy
    template <typename T2>
    //type( type<T2,tag,reference_object> && t ) : M_tmp( std::move( t.M_tmp ) ), M_value( std::move(std::forward<type<T2,tag,reference_object>>(t).M_value) ) {}
    //type( type<T2,tag,reference_object> && t ) : M_tmp( std::move( t.M_tmp ) ), M_value( std::move( static_cast<type<T2,tag,reference_object>&&>(t).M_value) ) {}
    //type( type<T2,tag,reference_object> && t ) : type( std::forward<type<T2,tag,reference_object>>(t).M_value ) {}
    constexpr type( type<T2,tag,reference_object> && t ) : type( std::move(t).M_value ) {}


    constexpr T & value()  { return M_value; }
    constexpr T const& value() const  { return M_value; }
    //T & value() && { return M_value; } // NOT ALWAYS GOOD!!
private :
    std::any M_tmp;
    T & M_value;

};




struct template_type_base : named_argument_base {};

template <typename tag>
struct template_type : template_type_base
{
    using tag_type = tag;

    template <typename T>
    using required_as_t = detail::infer_type_t<T,tag_type>;

    template <typename T>
    using apply_t = detail::infer_type_t<T,tag_type>;
};


} // detail




template <typename T, typename = void>
struct is_named_argument : std::false_type {};
// trick with decltype(sizeof(T)) to check if T is a complete type
template <typename T>
struct is_named_argument<T, std::void_t<decltype(sizeof(T))> > : std::is_base_of<detail::named_argument_base,T>::type {};

template <typename T>
constexpr bool is_named_argument_v = is_named_argument<T>::value;



template <typename TagType, typename /*Enable*/ = void>
struct named_argument_tag
{
    using type = TagType;
};
template <typename TagType>
struct named_argument_tag<TagType,std::void_t<std::enable_if_t< is_named_argument_v<TagType> , std::true_type/*bool*/>>>
{
    using type = typename TagType::tag_type;
};
template <typename TagType>
using named_argument_tag_t = typename named_argument_tag<TagType>::type;



template <typename NamedArgType, typename T>
struct named_argument
{
    using type = typename detail::infer_type_t<T, named_argument_tag_t<NamedArgType> >;

};
template <typename TagType>
struct named_argument<TagType,void>
{
    using type = detail::template_type<TagType>;

};
template <typename NamedArgType,typename T = void>
using named_argument_t = typename named_argument<NamedArgType,T>::type;



struct ArgumentIdentifierBase {};

template <typename T,typename /*Enable*/ = void>
struct ArgumentIdentifier : ArgumentIdentifierBase
{
    using identifier_type = T;
    constexpr ArgumentIdentifier() = default;
    ArgumentIdentifier(const ArgumentIdentifier&) = delete;
    ArgumentIdentifier& operator=(const ArgumentIdentifier&) = delete;

    template <typename V>
    constexpr typename T::template apply_t<V> operator=(V && v) const { return typename T::template apply_t<V>( std::forward<V>( v ) );}
};

template <typename T>
struct ArgumentIdentifier<T,std::void_t<std::enable_if_t< !is_named_argument_v<T> , std::true_type/*bool*/>>>  : ArgumentIdentifier<named_argument_t<T>,void>
{
    using ArgumentIdentifier<named_argument_t<T>,void>::operator=;
};


template <typename T>
constexpr ArgumentIdentifier<T> identifier;

template <typename T>
using identifier_t = ArgumentIdentifier<T>;

template <typename T>
constexpr bool is_argument_identifier_v = std::is_base_of_v<ArgumentIdentifierBase,std::decay_t<T>>;


//! make a new arg type with value \t
template <typename NamedArgType,typename T>
constexpr auto/*type<T,Tag>*/ make_argument(T&& t)
{
    return detail::infer_type_t<T,named_argument_tag_t<NamedArgType>>(std::forward<T>(t));
}

//! make a new arg type with value \t
template <typename ArgIdentifierType, typename T,  std::enable_if_t<is_argument_identifier_v<ArgIdentifierType> ,bool> = true >
constexpr auto make_argument(ArgIdentifierType &&, T&& t)
{
    return make_argument<typename std::decay_t<ArgIdentifierType>::identifier_type>( std::forward<T>(t) );
}



template <typename NamedArgType,typename ValueType>
struct DefaultArgument : named_argument_t<NamedArgType,ValueType>
{
    using super_type = named_argument_t<NamedArgType,ValueType>;
    DefaultArgument() = delete;
    template <typename U>
    constexpr DefaultArgument( U && v ) : named_argument_t<NamedArgType,ValueType>( std::forward<U>( v ) ) {}
    DefaultArgument( DefaultArgument const& ) = delete;
    DefaultArgument( DefaultArgument && ) = default;//delete;

    constexpr super_type const& to_named_argument() const & { return *this; }
    constexpr super_type & to_named_argument() & { return *this; }
    constexpr super_type && to_named_argument() && { return std::move(*this); }
};


template <typename NamedArgType,typename ValueType>
struct DefaultArgumentInvocable : named_argument_t<NamedArgType,ValueType>
{
    static_assert( std::is_invocable_v<ValueType>, "is_invocable_v should be true" );

    DefaultArgumentInvocable() = delete;
    template <typename U>
    constexpr DefaultArgumentInvocable( U && v ) : named_argument_t<NamedArgType,ValueType>( std::forward<U>( v ) ) {}
    DefaultArgumentInvocable( DefaultArgumentInvocable const& ) = delete;
    DefaultArgumentInvocable( DefaultArgumentInvocable && ) = default;//delete;

    constexpr auto to_named_argument() const & { return make_argument<NamedArgType>( this->value()() ); }
    constexpr auto to_named_argument() & { return make_argument<NamedArgType>( this->value()() ); }
    constexpr auto to_named_argument() && { return make_argument<NamedArgType>( std::move(*this).value()() ); }
};

template <typename NamedArgType,typename ValueType>
constexpr auto make_default_argument( ValueType && val )
{
    return DefaultArgument<NamedArgType,ValueType>{ std::forward<ValueType>( val ) };
}

template <typename ArgIdentifierType,typename ValueType, std::enable_if_t<is_argument_identifier_v<ArgIdentifierType> ,bool> = true >
constexpr auto make_default_argument( ArgIdentifierType &&, ValueType && val )
{
    return make_default_argument<typename std::decay_t<ArgIdentifierType>::identifier_type>( std::forward<ValueType>( val ) );
}


template <typename NamedArgType,typename ValueType>
constexpr auto make_default_argument_invocable( ValueType && val )
{
    return DefaultArgumentInvocable<NamedArgType,ValueType>{ std::forward<ValueType>( val ) };
}

template <typename ArgIdentifierType,typename ValueType, std::enable_if_t<is_argument_identifier_v<ArgIdentifierType> ,bool> = true >
constexpr auto make_default_argument_invocable( ArgIdentifierType &&, ValueType && val )
{
    return make_default_argument_invocable<typename std::decay_t<ArgIdentifierType>::identifier_type>( std::forward<ValueType>( val ) );
}


namespace detail
{

template<typename U, int Id, int MaxId, typename TupleType>
constexpr auto&& getImplBIS(TupleType && t)
{
    static_assert( Id< MaxId );
    if constexpr ( std::is_same_v<U, typename std::decay_t<std::tuple_element_t<Id, std::decay_t<TupleType> > >::tag_type > )
                     return std::get<Id>( std::forward<TupleType>(t) );//.value();
    else
        return getImplBIS<U,Id+1,MaxId,TupleType>( std::forward<TupleType>(t ));
}


#if 0
template<typename U, int Id, int MaxId, typename TupleType,typename DefaultArgsTupleType>
constexpr auto&& getImplBIS2(TupleType && t, DefaultArgsTupleType && defaultArgs)
{
    //static_assert( Id< MaxId );
    if constexpr ( Id>= MaxId )
                 {
                     return getImplBIS<U, 0, std::tuple_size<DefaultArgsTupleType>::value>( std::forward<DefaultArgsTupleType>( defaultArgs ) );
                 }
    else if constexpr ( std::is_same_v<U, typename std::tuple_element_t<Id, std::decay_t<TupleType> >::tag_type > )
                     return std::get<Id>( std::forward<TupleType>(t) );//.value();
    else
        return getImplBIS2<U,Id+1,MaxId,TupleType>( std::forward<TupleType>(t ), std::forward<DefaultArgsTupleType>( defaultArgs ) );
}

#endif

#if 0
template<typename U, int Id, int MaxId, typename DefaultArgsTupleType, typename RetTupleType>
constexpr auto /*decltype(auto)*/ newtuple( DefaultArgsTupleType && defaultArgs, RetTupleType && ret )
{
    if constexpr ( Id>= MaxId )
                     return std::forward<RetTupleType>( ret );
    else if constexpr ( U::template has_t< typename std::tuple_element_t<Id, std::decay_t<DefaultArgsTupleType> > >::value )
                          return newtuple<U,Id+1,MaxId>( std::forward<DefaultArgsTupleType>( defaultArgs ), std::forward<RetTupleType>( ret ) );
    else
        return newtuple<U,Id+1,MaxId>( std::forward<DefaultArgsTupleType>( defaultArgs ), std::tuple_cat( std::forward<RetTupleType>( ret ),
                                                                                                          std::/*make_tuple*/forward_as_tuple( std::get<Id>( std::forward<DefaultArgsTupleType>( defaultArgs ) ) )
                                       //std::make_tuple/*forward_as_tuple*/( std::get<Id>( defaultArgs  ) )
                                           //std::make_tuple/*forward_as_tuple*/( std::move(std::get<Id>( defaultArgs ) ) )
                                       )
                                       );
}
//auto newtuple = detail::subtuple<>( defaultArg... );
#else
template<typename U, int Id, int MaxId, typename DefaultArgsTupleType, typename ... RetTupleType>
constexpr decltype(auto) to_tuple_args_not_present( DefaultArgsTupleType && defaultArgs, RetTupleType && ... ret )
{
    if constexpr ( Id>= MaxId )
        return std::make_tuple( std::forward<RetTupleType>( ret )... );
    else if constexpr ( U::template has_t< typename std::tuple_element_t<Id, std::decay_t<DefaultArgsTupleType> > >::value )
        return to_tuple_args_not_present<U,Id+1,MaxId>( std::forward<DefaultArgsTupleType>( defaultArgs ), std::forward<RetTupleType>( ret )... );
    else
        return to_tuple_args_not_present<U,Id+1,MaxId>( std::forward<DefaultArgsTupleType>( defaultArgs ),
                                                        std::get<Id>( std::forward<DefaultArgsTupleType>( defaultArgs ) ).to_named_argument(), std::forward<RetTupleType>( ret )... );
}
#endif

template <typename ...T>
struct to_tuple : std::tuple<T...>
{
    using std::tuple<T...>::tuple;

    using super_type = std::tuple<T...>;
    // conversion operator
#if 0
    template <typename ... Ts>
    operator std::tuple<Ts...> () { return { std::get<Ts>(*this)...}; }
#else
    template <typename ... Ts>
    //operator std::tuple<Ts...> () { return std::tuple<Ts...>{ getImpl<typename Ts::tag_type,0, sizeof...(T) >( static_cast<super_type&&>( *this ) )...};}
    operator std::tuple<Ts...> () { return std::tuple<Ts...>{ getImplBIS<typename std::decay_t<Ts>::tag_type,0, sizeof...(T) >( static_cast<super_type&&>( *this ) )...}; }
#endif
};

template <typename ... T>
to_tuple(T...) -> to_tuple<T...>;

template<typename... Ts>
constexpr decltype(auto) to_tuple_from_tuple( std::tuple<Ts...> && theTuple)
{
    return std::apply
        (
            [](Ts &&... t) -> decltype(auto) { return to_tuple( std::forward<Ts>( t ) ... ); },
            std::forward<std::tuple<Ts...>>( theTuple )
         );
}


template <typename Tag,size_t index, class... Ts>
constexpr auto tuple_fold() {
    if constexpr( index >= sizeof...(Ts) )
                    return index;
    else if constexpr ( std::is_same_v<typename Tag::tag_type, typename std::decay_t<std::tuple_element_t<index, std::tuple<Ts...> > >::tag_type > )
                          return index;
    else
        return tuple_fold<Tag,index+1,Ts...>();

}
template <typename Tag,typename ... Ts>
constexpr auto tuple_index()
{
    return tuple_fold<Tag,0,Ts...>();
}

} // namespace detail


struct arguments_base {};

template <typename ... T>
struct arguments;

template <typename T>
constexpr bool is_arguments_v = std::is_base_of_v<arguments_base,std::decay_t<T>>;

template <typename ... T>
constexpr NA::arguments<T...>  make_arguments(T&& ... t);

template <typename ... Ts>
constexpr NA::arguments<Ts...>  make_arguments_from_tuple( std::tuple<Ts...> &&  t )
{
    return std::apply
        (
            [](Ts &&... t) -> decltype(auto) { return make_arguments( std::forward<Ts>( t ) ... ); },
            //[](auto &&... t) -> decltype(auto) { return make_arguments( std::forward<decltype(auto)>( t ) ... ); },
            std::forward<std::tuple<Ts...>>( t )
         );

}


template <typename ... T>
struct arguments : arguments_base
{
    // using tuple_type = std::tuple<T...>;

    template <typename Tag>
    using has_t = std::disjunction<std::is_same<typename std::decay_t<T>::tag_type, typename std::decay_t<Tag>::tag_type>...>;

    template <typename ... U, typename  = typename std::enable_if_t< sizeof...(U) >= 1 > >
    constexpr arguments(U&& ... u) : values(detail::to_tuple(std::forward<U>(u)...)) {}

    template <typename U, std::enable_if_t< !is_arguments_v<U>, bool > = true >
    constexpr arguments(U&&  u) : values(std::forward<U>(u)) {}

    template <typename U, std::enable_if_t< is_arguments_v<U>, bool > = true >
    constexpr arguments(U&& u) : values( to_tuple_from_tuple( std::forward<U>(u).values ) ) {}

    //arguments() = delete;
    template< typename T2 = std::tuple<T...>,
              typename  = typename std::enable_if_t< std::tuple_size<T2>::value == 0 > >
    constexpr arguments() {}

    arguments( arguments const& ) = delete;
    //arguments( arguments && ) = delete;


    template <typename NamedArgType>
    constexpr auto const& get() const { return this->getArgument<NamedArgType>().value(); }

    template <typename NamedArgType>
    constexpr auto & get() { return this->getArgument<NamedArgType>().value(); }

    template <typename ArgIdentifierType,  std::enable_if_t<is_argument_identifier_v<ArgIdentifierType> ,bool> = true >
    constexpr auto const& get( ArgIdentifierType && t ) const { return this->get<typename std::decay_t<ArgIdentifierType>::identifier_type>(); }

    template <typename ArgIdentifierType,  std::enable_if_t<is_argument_identifier_v<ArgIdentifierType> ,bool> = true >
    constexpr auto & get( ArgIdentifierType && t ) { return this->get<typename std::decay_t<ArgIdentifierType>::identifier_type>(); }

    // check constraint of return type
    template <template <typename RT> class ReturnConstraintType, typename ArgIdentifierType,  std::enable_if_t<is_argument_identifier_v<ArgIdentifierType> ,bool> = true >
    constexpr auto const& get( ArgIdentifierType && t ) const
        {
            static_assert( ReturnConstraintType< std::decay_t<decltype(this->get<typename std::decay_t<ArgIdentifierType>::identifier_type>())> >::value,
                           "argument type not statisfy the constraint asked" );
            return this->get<typename std::decay_t<ArgIdentifierType>::identifier_type>();
        }

    template <template <typename RT> class ReturnConstraintType, typename ArgIdentifierType,  std::enable_if_t<is_argument_identifier_v<ArgIdentifierType> ,bool> = true >
    constexpr auto & get( ArgIdentifierType && t )
        {
            static_assert( ReturnConstraintType< std::decay_t<decltype(this->get<typename std::decay_t<ArgIdentifierType>::identifier_type>())> >::value,
                           "argument type not statisfy the constraint asked" );
            return this->get<typename std::decay_t<ArgIdentifierType>::identifier_type>();
        }


    //---------------------------------
    // get_else
    //---------------------------------
    // const
    template <typename NamedArgType,typename DefaultType, std::enable_if_t< has_t<NamedArgType>::value || std::is_lvalue_reference_v<DefaultType&&>, bool > = true >
    constexpr decltype(auto) /*auto &&*/ get_else( DefaultType && defaultValue ) const
        {
            if constexpr ( has_t<NamedArgType>::value )
                return this->get<NamedArgType>();
            else
                return std::forward<DefaultType>( defaultValue );
        }

    template <typename NamedArgType,typename DefaultType, std::enable_if_t< !has_t<NamedArgType>::value && std::is_rvalue_reference_v<DefaultType&&>, bool > = true >
    constexpr auto get_else( DefaultType && defaultValue ) const
        {
            //return defaultValue;
            return std::forward<DefaultType>( defaultValue );
        }

    template <typename ArgIdentifierType, typename DefaultType>
    constexpr decltype(auto) get_else( ArgIdentifierType && t, DefaultType && defaultValue ) const
        {
            return this->get_else<typename std::decay_t<ArgIdentifierType>::identifier_type>( std::forward<DefaultType>( defaultValue ) );
        }
    // non const
    template <typename NamedArgType,typename DefaultType, std::enable_if_t< has_t<NamedArgType>::value || std::is_lvalue_reference_v<DefaultType&&>, bool > = true >
    constexpr decltype(auto) /*auto &&*/ get_else( DefaultType && defaultValue )
        {
            if constexpr ( has_t<NamedArgType>::value )
                return this->get<NamedArgType>();
            else
                return std::forward<DefaultType>( defaultValue );
        }
    template <typename ArgIdentifierType, typename DefaultType>
    constexpr decltype(auto) get_else( ArgIdentifierType && t, DefaultType && defaultValue )
        {
            return this->get_else<typename std::decay_t<ArgIdentifierType>::identifier_type>( std::forward<DefaultType>( defaultValue ) );
        }
    // constraint
    template <template <typename RT> class ReturnConstraintType, typename ArgIdentifierType, typename DefaultType >
    constexpr decltype(auto) get_else( ArgIdentifierType && t, DefaultType && defaultValue ) const
        {
            decltype(auto) res = this->get_else<typename std::decay_t<ArgIdentifierType>::identifier_type>( std::forward<DefaultType>( defaultValue ) );
            static_assert( ReturnConstraintType< std::decay_t<decltype(res)> >::value,
                           "argument type not statisfy the constraint asked" );
            return res;
        }

    //---------------------------------
    // get_else_invocable
    //---------------------------------
    // const
    template <typename NamedArgType,typename DefaultType, std::enable_if_t< has_t<NamedArgType>::value || std::is_lvalue_reference_v<  std::invoke_result_t< DefaultType> >, bool > = true >
    constexpr decltype(auto) /*auto &&*/ get_else_invocable( DefaultType && defaultInvocableValue ) const
        {
            if constexpr ( has_t<NamedArgType>::value )
                return this->get<NamedArgType>();
            else
            {
                static_assert( std::is_invocable_v<DefaultType>, "is_invocable_v should be true" );
                return std::forward<DefaultType>( defaultInvocableValue )();
            }
        }
    template <typename NamedArgType,typename DefaultType,
              std::enable_if_t< !has_t<NamedArgType>::value && ( !std::is_reference_v< std::invoke_result_t< DefaultType> > || std::is_rvalue_reference_v< std::invoke_result_t< DefaultType> > ) , bool > = true >
    constexpr auto get_else_invocable( DefaultType && defaultInvocableValue ) const
        {
            static_assert( std::is_invocable_v<DefaultType>, "is_invocable_v should be true" );
            return std::forward<DefaultType>( defaultInvocableValue )();
        }
    template <typename ArgIdentifierType,typename DefaultType>
    constexpr decltype(auto) get_else_invocable( ArgIdentifierType && t, DefaultType && defaultInvocableValue ) const
        {
            return this->get_else_invocable<typename std::decay_t<ArgIdentifierType>::identifier_type>( std::forward<DefaultType>( defaultInvocableValue ) );
        }
    // non const
    template <typename NamedArgType,typename DefaultType, std::enable_if_t< has_t<NamedArgType>::value || std::is_lvalue_reference_v<  std::invoke_result_t< DefaultType> >, bool > = true >
    constexpr decltype(auto) /*auto &&*/ get_else_invocable( DefaultType && defaultInvocableValue )
        {
            if constexpr ( has_t<NamedArgType>::value )
                return this->get<NamedArgType>();
            else
            {
                static_assert( std::is_invocable_v<DefaultType>, "is_invocable_v should be true" );
                return std::forward<DefaultType>( defaultInvocableValue )();
            }
        }
    template <typename ArgIdentifierType,typename DefaultType>
    constexpr decltype(auto) get_else_invocable( ArgIdentifierType && t, DefaultType && defaultInvocableValue )
        {
            return this->get_else_invocable<typename std::decay_t<ArgIdentifierType>::identifier_type>( std::forward<DefaultType>( defaultInvocableValue ) );
        }

    // with constraint
    template <template <typename RT> class ReturnConstraintType, typename ArgIdentifierType,typename DefaultType>
    constexpr decltype(auto) get_else_invocable( ArgIdentifierType && t, DefaultType && defaultInvocableValue ) const
        {
            decltype(auto) res = this->get_else_invocable<typename std::decay_t<ArgIdentifierType>::identifier_type>( std::forward<DefaultType>( defaultInvocableValue ) );
            static_assert( ReturnConstraintType< std::decay_t<decltype(res)> >::value,
                           "argument type not statisfy the constraint asked" );
            return res;
        }


    //! return all argument
    constexpr auto get_all() const { return std::forward_as_tuple( this->get<std::decay_t<T>>() ... ); }
    //! return all argument
    constexpr auto get_all() { return std::forward_as_tuple( this->get<std::decay_t<T>>() ... ); }


    template <typename NamedArgType>
    constexpr auto & getArgument() &
        {
            //return std::get< detail::tuple_index<NamedArgType,T...>() >( this->values );
            if constexpr ( has_t<NamedArgType>::value ) // non template type
                return std::get< detail::tuple_index<NamedArgType,T...>() >( this->values );
            else
                return detail::getImplBIS<typename NamedArgType::tag_type,0,sizeof...(T)>( this->values );
        }
    template <typename NamedArgType>
    constexpr auto const& getArgument() const &
    {
        if constexpr ( false ) // non template type
            return std::get<NamedArgType>( this->values ).value;
        else
            return detail::getImplBIS<typename NamedArgType::tag_type,0,sizeof...(T)>( this->values );
    }
    template <typename NamedArgType>
    constexpr auto && getArgument() &&
        {

            if constexpr ( false ) // non template type
                return std::get<NamedArgType>( std::move( this->values ) );
            else
                return detail::getImplBIS<typename NamedArgType::tag_type,0,sizeof...(T)>( std::move( this->values ) );
        }


#if 0
    template <typename NamedArgType, typename DefaultArgsTupleType>
    constexpr auto && getOptionalArgument( DefaultArgsTupleType && defaultArgs ) //&&
    {
        return detail::getImplBIS2<typename NamedArgType::tag_type,0,sizeof...(T)>( std::move( this->values ), std::forward<DefaultArgsTupleType>( defaultArgs )  );
    }
#endif

    template <typename ... DefaultArgType >
    constexpr decltype(auto) add_default_arguments( DefaultArgType && ... defaultArg ) &&
        {
            //static constexpr bool __matches[sizeof...(_Args)] = {is_same<_T1, _Args>::value...};
            using self_type = arguments<T...>;
            auto tupleDefaultArgsUsed = detail::to_tuple_args_not_present<self_type,0,sizeof...(DefaultArgType)>( std::make_tuple( std::forward<DefaultArgType>(defaultArg)... ) );

            if constexpr ( std::tuple_size<std::decay_t<decltype(tupleDefaultArgsUsed)>>::value == 0 )
                return std::move( *this );
            else
                return make_arguments_from_tuple( std::tuple_cat( std::move( this->values), std::move( tupleDefaultArgsUsed ) ) );
        }

    template <typename GivenArgsType, typename ... DefaultArgType >
    static constexpr auto create( GivenArgsType && givenArgs, DefaultArgType && ... defaultArg )
        {
            //return createImpl( std::forward<GivenArgsType>( givenArgs ), std::make_tuple( std::forward<DefaultArgType>( defaultArg ) ... ) );

            return arguments<T...>{ std::forward<GivenArgsType>( givenArgs).add_default_arguments( std::forward<DefaultArgType>( defaultArg )  ... ) };
        }

#if 0
private :
    template <typename TheType, typename DefaultArgsTupleType >
    static constexpr auto createImpl( TheType && a, DefaultArgsTupleType &&  defaultParams )
        {
            //return args<T...>{ std::forward<TheType>( a ).template getArgument<T>() ... };
#if 1
            return arguments<T...>{ std::forward<TheType>( a ).template getOptionalArgument<T>( std::forward<DefaultArgsTupleType>( defaultParams ) ) ... };
#endif
            //return arguments<T...>{ std::forward<TheType>( a ).add_default_arguments( std::forward<DefaultArgsTupleType>( defaultParams )  ... ) };
        }
#endif

private :

    template <typename ... T2> friend struct arguments;

    std::tuple<T...> values;
};




template <typename ... T>
constexpr NA::arguments<T...>  make_arguments(T&& ... t)
{
    return NA::arguments<T...>{ std::forward<T>(t)... };
}

template <typename U,typename ... T>
constexpr U  make_arguments(T&& ... t)
{
    static_assert( is_arguments_v<U> );
    return U{ std::forward<T>(t)... };
}


template <typename Tag,typename U>
constexpr bool has_v = U::template has_t<Tag>::value;





namespace constraint
{

template <typename To>
struct is_convertible
{
    template <typename From>
    struct apply : std::is_convertible<From,To> {};
};

} // constraint





} //  namespace NA

#endif