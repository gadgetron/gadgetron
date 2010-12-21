#include <ace/Log_Msg.h>
#include <ace/Service_Config.h>
#include <ace/Reactor.h>


int ACE_TMAIN(int argc, ACE_TCHAR *argv[])
{
  ACE_TRACE(( ACE_TEXT("main") ));

  ACE_LOG_MSG->priority_mask( LM_INFO | LM_NOTICE | LM_ERROR| LM_DEBUG,
			     ACE_Log_Msg::PROCESS);

  ACE_DEBUG(( LM_DEBUG, ACE_TEXT("%IConfiguring services\n") ));

  //ACE_Service_Config::open(argc, argv);

  ACE_Reactor::instance()->run_reactor_event_loop ();

  return 0;
}
