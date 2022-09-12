//
// Created by dchansen on 4/3/19.
//

#pragma once
#include "MPMCChannel.h"
#include <boost/hana.hpp>
#include <future>

namespace Gadgetron::Core {


    class ThreadPool {
    private:
        class Work {
        public:
            virtual void execute() = 0;
            virtual ~Work()        = default;
        };

        template <class F, class... ARGS> class Storage{
        public:
            explicit Storage(F&& f, ARGS&&... args)
                : f{ std::forward<F>(f) }, args{ std::forward<ARGS>(args)... } {}
        protected:
            F f;
            boost::hana::tuple<std::decay_t<ARGS>...> args;
        public:
            using R = decltype(boost::hana::unpack(std::move(args),f));
        protected:
            std::promise<R> promise;


        public:
            std::future<R> get_future() {
                return promise.get_future();
            }
        };
        template<class F, bool IS_VOID, class... ARGS> class ConcreteWorkImpl{};

        template<class F, class... ARGS> class ConcreteWorkImpl<F,false,ARGS...> : public Work, public Storage<F,ARGS...> {
        public:
            using Storage<F,ARGS...>::Storage;

            void execute() override {
                try {
                    auto val = boost::hana::unpack(std::move(this->args),this->f);
                    this->promise.set_value(std::move(val));
                } catch (...){
                    this->promise.set_exception(std::current_exception());
                }
            }
        };
        template<class F, class... ARGS> class ConcreteWorkImpl<F,true,ARGS...> : public Work, public Storage<F,ARGS...> {
        public:
            using Storage<F,ARGS...>::Storage;

            void execute() override {
                try {
                    boost::hana::unpack(std::move(this->args),this->f);
                    this->promise.set_value();
                } catch (...){
                    this->promise.set_exception(std::current_exception());
                }
            }
        };

    template<class F, class... ARGS> class ConcreteWork : public ConcreteWorkImpl<F, std::is_same<typename Storage<F,ARGS...>::R,void>::value,ARGS...> {
        public:
            using ConcreteWorkImpl<F, std::is_same<typename Storage<F,ARGS...>::R,void>::value,ARGS...>::ConcreteWorkImpl;
        };

    public:
        explicit ThreadPool(unsigned int workers) {
            for (auto i = 0u; i < workers; i++) {
                threads.emplace_back([this]() {
                    try {
                        while (true)
                            this->work_queue.pop()->execute();
                    } catch (const ChannelClosed&) {
                    }
                });
            }
        }

        template <class F, class... ARGS> auto async(F&& f, ARGS&&... args) {
            auto work = std::make_unique<ConcreteWork<F, ARGS...>>(std::forward<F>(f), std::forward<ARGS>(args)...);
            auto future_result = work->get_future();
            work_queue.push(std::move(work));
            return std::move(future_result);
        }

        void join(){
            work_queue.close();
            for (auto& thread : threads){
                thread.join();
            }
        }

    private:
        MPMCChannel<std::unique_ptr<Work>> work_queue;
        std::vector<std::thread> threads;
    };

}
