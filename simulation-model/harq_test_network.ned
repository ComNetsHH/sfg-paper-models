import inet.networklayer.configurator.ipv4.IPv4NetworkConfigurator;
import inet.networklayer.ipv4.RoutingTableRecorder;
import inet.node.ethernet.Eth10G;
import inet.node.inet.Router;
import inet.node.inet.StandardHost;
import lte.corenetwork.binder.LteBinder;
import lte.corenetwork.nodes.Ue;
import lte.corenetwork.nodes.eNodeB;
import lte.epc.PgwStandardSimplified;
import lte.world.radio.LteChannelControl;


network harq_test_network
{
    @display("bgb=1144,754");
    submodules:
        channelControl: LteChannelControl {
            @display("p=408,227");
        }
        routingTableRecorder: RoutingTableRecorder {
            @display("p=603,228");
        }
        binder: LteBinder {
            @display("p=505,228");
        }
        pgw: PgwStandardSimplified {
            nodeType = "PGW";
            @display("p=13,107");
        }
        eNodeB: eNodeB {
            @display("p=13,20");
        }
        ue: Ue {
            @display("p=185,20");
        }
        router: Router {
            @display("p=13,185");
        }
        server: StandardHost {
            @display("p=13,275");
        }
        configurator: IPv4NetworkConfigurator {
            @display("p=703,234");
        }
    connections:
        pgw.pppg++ <--> Eth10G <--> eNodeB.ppp;
        pgw.filterGate <--> Eth10G <--> router.pppg++;
        router.pppg++ <--> Eth10G <--> server.pppg++;
}
